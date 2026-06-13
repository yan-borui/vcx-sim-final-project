#include "Labs/FinalProject/FreeSurfaceSeparationSimulator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

namespace VCX::Labs::Final {
    namespace {
        struct FaceNeighbor {
            glm::ivec3 CellOffset;
            glm::ivec3 FaceOffset;
            int        Direction;
            float      DivergenceSign;
        };

        constexpr std::array<FaceNeighbor, 6> FaceNeighbors {
            FaceNeighbor { { -1, 0, 0 }, { 0, 0, 0 }, 0, -1.0f },
            FaceNeighbor {  { 1, 0, 0 }, { 1, 0, 0 }, 0,  1.0f },
            FaceNeighbor { { 0, -1, 0 }, { 0, 0, 0 }, 1, -1.0f },
            FaceNeighbor {  { 0, 1, 0 }, { 0, 1, 0 }, 1,  1.0f },
            FaceNeighbor { { 0, 0, -1 }, { 0, 0, 0 }, 2, -1.0f },
            FaceNeighbor {  { 0, 0, 1 }, { 0, 0, 1 }, 2,  1.0f },
        };

        // The standalone case uses grid-aligned tank walls, exactly half a cell
        // from the adjacent pressure sample.
        constexpr float TankWallGhostScale = 2.0f;
    }

    void FreeSurfaceSeparationSimulator::setupScene(int res) {
        Simulator::setupScene(res);

        m_body            = nullptr;
        gravity           = { 0.0f, -9.81f, 0.0f };
        m_fRatio          = 0.97f;
        numSubSteps       = 2;
        numParticleIters  = 3;
        numPressureIters  = 120;
        separateParticles = true;
        compensateDrift   = false;

        setupSplashParticles();
        wallSeparationCandidates          = 0;
        separatingWallCells               = 0;
        wallSeparationActiveSetIterations = 0;
        separatingWallFaces               = 0;
        wallContactParticles              = 0;
        minimumContactPressure            = 0.0f;
        pressureResidual                  = 0.0f;
        wallSeparationKktResidual         = 0.0f;
        averageLeftWallDistance           = 0.0f;
        updateWallContactDiagnostics();
    }

    void FreeSurfaceSeparationSimulator::setupSplashParticles() {
        // A thin sheet makes the wall attachment visible from the default
        // front camera, matching the 2D comparison in Figure 6 of the paper.
        glm::vec3 const center { 0.12f, 0.16f, 0.0f };
        glm::vec3 const radius { 0.18f, 0.18f, 0.08f };
        float const     spacing = std::max(2.0f * m_particleRadius, 0.42f * m_h);

        std::vector<glm::vec3> positions;
        for (float x = center.x - radius.x; x <= center.x + radius.x; x += spacing) {
            int const xLayer = int(std::round((x - center.x + radius.x) / spacing));
            for (float y = center.y - radius.y; y <= center.y + radius.y; y += spacing) {
                int const   yLayer  = int(std::round((y - center.y + radius.y) / spacing));
                float const zOffset = ((xLayer + yLayer) & 1) ? 0.5f * spacing : 0.0f;
                for (float z = center.z - radius.z + zOffset; z <= center.z + radius.z; z += spacing) {
                    glm::vec3 const normalized = (glm::vec3(x, y, z) - center) / radius;
                    if (glm::dot(normalized, normalized) <= 1.0f)
                        positions.emplace_back(x, y, z);
                }
            }
        }

        m_iNumSpheres = int(positions.size());
        m_particlePos = std::move(positions);
        m_particleVel.assign(m_iNumSpheres, glm::vec3(-3.0f, 0.7f, 0.0f));
        m_particleColor.assign(m_iNumSpheres, glm::vec3(1.0f));
        m_hashtable.assign(m_iNumSpheres, 0);
        m_hashtableindex.assign(m_iNumCells + 1, 0);
        m_particleRestDensity = 0.0f;
    }

    void FreeSurfaceSeparationSimulator::SimulateTimestep(float dt) {
        if (dt <= 0.0f)
            return;

        float const substep = dt / float(std::max(numSubSteps, 1));
        for (int step = 0; step < std::max(numSubSteps, 1); ++step) {
            integrateParticles(substep);
            constrainParticlesToTank();
            if (separateParticles)
                pushParticlesApart(numParticleIters);
            constrainParticlesToTank();
            transferVelocities(true, m_fRatio);
            updateParticleDensity();
            solveIncompressibility(numPressureIters, substep, overRelaxation, compensateDrift);
            transferVelocities(false, m_fRatio);
        }
        updateParticleColors();
        updateWallContactDiagnostics();
    }

    void FreeSurfaceSeparationSimulator::constrainParticlesToTank() {
        float const minBound = -0.5f + m_h + m_particleRadius;
        float const maxBound = 0.5f - m_h - m_particleRadius;
        for (glm::vec3 & position : m_particlePos)
            position = glm::clamp(position, glm::vec3(minBound), glm::vec3(maxBound));
    }

    bool FreeSurfaceSeparationSimulator::isValidCell(glm::ivec3 const & cell) const {
        return cell.x >= 0 && cell.x < m_iCellX
            && cell.y >= 0 && cell.y < m_iCellY
            && cell.z >= 0 && cell.z < m_iCellZ;
    }

    int FreeSurfaceSeparationSimulator::gridOffset(glm::ivec3 const & cell) const {
        return cell.x + m_iCellX * (cell.y + m_iCellY * cell.z);
    }

    bool FreeSurfaceSeparationSimulator::isSolidCell(glm::ivec3 const & cell) const {
        return ! isValidCell(cell) || m_s[gridOffset(cell)] <= 0.0f;
    }

    bool FreeSurfaceSeparationSimulator::isWallSeparationCandidate(glm::ivec3 const & cell) const {
        for (FaceNeighbor const & side : FaceNeighbors) {
            glm::ivec3 const neighbor = cell + side.CellOffset;
            if (isSolidCell(neighbor))
                return true;
        }
        return false;
    }

    void FreeSurfaceSeparationSimulator::solveIncompressibility(
        int   numIters,
        float dt,
        float overRelaxation,
        bool  compensateDrift) {
        (void) dt;
        (void) overRelaxation;
        (void) compensateDrift;

        std::vector<float> const liquidLevelSet = buildParticleLevelSet();
        std::vector<int> cellToRow(m_iNumCells, -1);
        std::vector<int> rowToCell;
        rowToCell.reserve(m_iNumCells / 3);
        for (int idx = 0; idx < m_iNumCells; ++idx) {
            if (m_type[idx] == FLUID_CELL) {
                cellToRow[idx] = int(rowToCell.size());
                rowToCell.push_back(idx);
            }
        }

        std::fill(m_p.begin(), m_p.end(), 0.0f);
        wallSeparationCandidates          = 0;
        separatingWallCells               = 0;
        wallSeparationActiveSetIterations = 0;
        separatingWallFaces               = 0;
        minimumContactPressure            = 0.0f;
        pressureResidual                  = 0.0f;
        wallSeparationKktResidual         = 0.0f;
        if (rowToCell.empty())
            return;

        auto decodeCell = [&](int idx) {
            return glm::ivec3(
                idx % m_iCellX,
                (idx / m_iCellX) % m_iCellY,
                idx / (m_iCellX * m_iCellY));
        };

        struct WallFace {
            int   Row;
            int   FaceIndex;
            int   Direction;
            float DivergenceSign;
            float WallVelocity;
            bool  Candidate;
            bool  Contact;
        };

        std::vector<WallFace> wallFaces;
        std::vector<int>      wallFaceForSide(rowToCell.size() * FaceNeighbors.size(), -1);
        for (int row = 0; row < int(rowToCell.size()); ++row) {
            glm::ivec3 const cell          = decodeCell(rowToCell[row]);
            bool const       candidateCell = isWallSeparationCandidate(cell);
            if (candidateCell)
                ++wallSeparationCandidates;

            for (int sideIndex = 0; sideIndex < int(FaceNeighbors.size()); ++sideIndex) {
                FaceNeighbor const & side     = FaceNeighbors[sideIndex];
                glm::ivec3 const     neighbor = cell + side.CellOffset;
                if (! isSolidCell(neighbor))
                    continue;

                glm::ivec3 const face = cell + side.FaceOffset;
                wallFaceForSide[row * FaceNeighbors.size() + sideIndex] =
                    int(wallFaces.size());
                wallFaces.push_back(WallFace {
                    .Row            = row,
                    .FaceIndex      = gridOffset(face),
                    .Direction      = side.Direction,
                    .DivergenceSign = side.DivergenceSign,
                    .WallVelocity   = 0.0f,
                    .Candidate      = enableWallSeparation && candidateCell,
                    .Contact        = true,
                });
            }
        }

        std::vector<glm::vec3> const intermediateVelocity = m_vel;
        for (WallFace & wallFace : wallFaces) {
            if (! wallFace.Candidate)
                continue;

            float const unconstrainedSeparationSpeed =
                -wallFace.DivergenceSign
                * (intermediateVelocity[wallFace.FaceIndex][wallFace.Direction]
                   - wallFace.WallVelocity);
            if (unconstrainedSeparationSpeed > 1e-5f)
                wallFace.Contact = false;
        }

        Eigen::VectorXd              pressure             = Eigen::VectorXd::Zero(int(rowToCell.size()));
        auto const contactPressure = [&](WallFace const & wallFace) {
            float const velocityDelta =
                wallFace.WallVelocity
                - intermediateVelocity[wallFace.FaceIndex][wallFace.Direction];
            return float(pressure[wallFace.Row])
                - velocityDelta
                    / (wallFace.DivergenceSign * TankWallGhostScale);
        };
        auto                         updateMinimumContactPressure = [&]() {
            bool  hasContact = false;
            float minimum    = std::numeric_limits<float>::infinity();
            for (WallFace const & wallFace : wallFaces) {
                if (! wallFace.Candidate || ! wallFace.Contact)
                    continue;
                hasContact = true;
                minimum    = std::min(minimum, contactPressure(wallFace));
            }
            minimumContactPressure = hasContact ? minimum : 0.0f;
        };
        auto                         solveWithContacts    = [&]() {
            int const                           matrixSize = int(rowToCell.size());
            Eigen::SparseMatrix<double>         matrix(matrixSize, matrixSize);
            std::vector<Eigen::Triplet<double>> triplets;
            triplets.reserve(matrixSize * 7);
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(matrixSize);

            std::vector<char> pinnedRows(matrixSize, 0);
            std::vector<char> visitedRows(matrixSize, 0);
            for (int seed = 0; seed < matrixSize; ++seed) {
                if (visitedRows[seed])
                    continue;

                bool             hasDirichletBoundary = false;
                std::vector<int> component;
                std::vector<int> pending { seed };
                visitedRows[seed] = 1;
                while (! pending.empty()) {
                    int const row = pending.back();
                    pending.pop_back();
                    component.push_back(row);

                    glm::ivec3 const cell = decodeCell(rowToCell[row]);
                    for (int sideIndex = 0; sideIndex < int(FaceNeighbors.size()); ++sideIndex) {
                        FaceNeighbor const & side     = FaceNeighbors[sideIndex];
                        glm::ivec3 const     neighbor = cell + side.CellOffset;
                        int const            wallFaceIndex =
                            wallFaceForSide[row * FaceNeighbors.size() + sideIndex];
                        if (wallFaceIndex >= 0) {
                            WallFace const & wallFace = wallFaces[wallFaceIndex];
                            if (wallFace.Candidate && ! wallFace.Contact)
                                hasDirichletBoundary = true;
                            continue;
                        }

                        int const neighborRow = cellToRow[gridOffset(neighbor)];
                        if (neighborRow >= 0) {
                            if (! visitedRows[neighborRow]) {
                                visitedRows[neighborRow] = 1;
                                pending.push_back(neighborRow);
                            }
                        } else if (m_type[gridOffset(neighbor)] == EMPTY_CELL) {
                            hasDirichletBoundary = true;
                        }
                    }
                }

                if (! hasDirichletBoundary)
                    pinnedRows[component.front()] = 1;
            }

            for (int row = 0; row < matrixSize; ++row) {
                if (pinnedRows[row]) {
                    triplets.emplace_back(row, row, 1.0);
                    continue;
                }

                glm::ivec3 const cell       = decodeCell(rowToCell[row]);
                double           diagonal   = 0.0;
                double           divergence = 0.0;
                for (int sideIndex = 0; sideIndex < int(FaceNeighbors.size()); ++sideIndex) {
                    FaceNeighbor const & side     = FaceNeighbors[sideIndex];
                    glm::ivec3 const     face     = cell + side.FaceOffset;
                    glm::ivec3 const     neighbor = cell + side.CellOffset;
                    int const            faceIdx  = gridOffset(face);
                    int const            wallFaceIndex =
                        wallFaceForSide[row * FaceNeighbors.size() + sideIndex];

                    if (wallFaceIndex >= 0) {
                        WallFace const & wallFace = wallFaces[wallFaceIndex];
                        float const      velocity =
                            wallFace.Candidate && ! wallFace.Contact
                            ? intermediateVelocity[faceIdx][side.Direction]
                            : wallFace.WallVelocity;
                        divergence += side.DivergenceSign * double(velocity);
                        if (wallFace.Candidate && ! wallFace.Contact)
                            diagonal += TankWallGhostScale;
                        continue;
                    }

                    divergence += side.DivergenceSign
                        * double(intermediateVelocity[faceIdx][side.Direction]);
                    int const neighborIdx = gridOffset(neighbor);
                    int const neighborRow = cellToRow[neighborIdx];
                    if (neighborRow >= 0) {
                        diagonal += 1.0;
                        if (! pinnedRows[neighborRow])
                            triplets.emplace_back(row, neighborRow, -1.0);
                    } else if (m_type[neighborIdx] == EMPTY_CELL) {
                        diagonal += ghostFluidPressureScale(
                            liquidLevelSet,
                            rowToCell[row],
                            neighborIdx);
                    } else {
                        diagonal += 1.0;
                    }
                }

                if (diagonal <= 1e-8) {
                    diagonal   = 1.0;
                    divergence = 0.0;
                }
                triplets.emplace_back(row, row, diagonal);
                rhs[row] = -divergence;
            }

            matrix.setFromTriplets(triplets.begin(), triplets.end());
            Eigen::ConjugateGradient<
                Eigen::SparseMatrix<double>,
                Eigen::Lower | Eigen::Upper,
                Eigen::IncompleteCholesky<double>>
                solver;
            solver.setMaxIterations(std::max(numIters, 1));
            solver.setTolerance(1e-5);
            solver.compute(matrix);
            if (solver.info() == Eigen::Success)
                pressure = solver.solve(rhs);
            if (solver.info() == Eigen::Success && pressure.allFinite()) {
                pressureResidual = float(solver.error());
            } else {
                pressure.setZero();
                pressureResidual = std::numeric_limits<float>::infinity();
            }

            m_vel = intermediateVelocity;
            for (WallFace const & wallFace : wallFaces) {
                if (! wallFace.Candidate || wallFace.Contact)
                    m_vel[wallFace.FaceIndex][wallFace.Direction] =
                        wallFace.WallVelocity;
            }

            for (int row = 0; row < matrixSize; ++row) {
                glm::ivec3 const cell         = decodeCell(rowToCell[row]);
                float const      cellPressure = float(pressure[row]);
                for (int sideIndex = 0; sideIndex < int(FaceNeighbors.size()); ++sideIndex) {
                    FaceNeighbor const & side    = FaceNeighbors[sideIndex];
                    glm::ivec3 const     face    = cell + side.FaceOffset;
                    glm::ivec3 const     neighbor = cell + side.CellOffset;
                    int const            faceIdx = gridOffset(face);
                    float                pressureScale = 1.0f;
                    int const            wallFaceIndex =
                        wallFaceForSide[row * FaceNeighbors.size() + sideIndex];
                    if (wallFaceIndex >= 0) {
                        WallFace const & wallFace = wallFaces[wallFaceIndex];
                        if (! wallFace.Candidate || wallFace.Contact)
                            continue;
                        pressureScale = TankWallGhostScale;
                    } else if (m_type[gridOffset(neighbor)] == EMPTY_CELL) {
                        pressureScale = ghostFluidPressureScale(
                            liquidLevelSet,
                            rowToCell[row],
                            gridOffset(neighbor));
                    }
                    m_vel[faceIdx][side.Direction] +=
                        side.DivergenceSign * pressureScale * cellPressure;
                }
            }
        };

        bool      changedActiveSet       = false;
        int const maxActiveSetIterations = enableWallSeparation
            ? std::min<int>(256, 2 * int(wallFaces.size()) + 8)
            : 1;
        for (int iter = 0; iter < maxActiveSetIterations; ++iter) {
            solveWithContacts();
            wallSeparationActiveSetIterations = iter + 1;
            updateMinimumContactPressure();

            changedActiveSet = false;
            if (! enableWallSeparation)
                break;

            bool releasedSuctionContact = false;
            for (WallFace & wallFace : wallFaces) {
                if (! wallFace.Candidate || ! wallFace.Contact)
                    continue;

                float const boundaryPressure = contactPressure(wallFace);
                if (boundaryPressure < -1e-5f) {
                    wallFace.Contact       = false;
                    releasedSuctionContact = true;
                    changedActiveSet       = true;
                }
            }
            if (releasedSuctionContact)
                continue;

            for (WallFace & wallFace : wallFaces) {
                if (! wallFace.Candidate || wallFace.Contact)
                    continue;

                float const separationSpeed =
                    -wallFace.DivergenceSign
                    * (m_vel[wallFace.FaceIndex][wallFace.Direction]
                       - wallFace.WallVelocity);
                if (separationSpeed < -1e-5f) {
                    wallFace.Contact = true;
                    changedActiveSet = true;
                }
            }
            if (! changedActiveSet)
                break;
        }

        if (changedActiveSet) {
            solveWithContacts();
            ++wallSeparationActiveSetIterations;
            updateMinimumContactPressure();
        }

        std::vector<char> separatingCell(rowToCell.size(), 0);
        separatingWallFaces = 0;
        wallSeparationKktResidual = 0.0f;
        for (WallFace const & wallFace : wallFaces) {
            if (! wallFace.Candidate)
                continue;

            float const separationSpeed =
                -wallFace.DivergenceSign
                * (m_vel[wallFace.FaceIndex][wallFace.Direction]
                   - wallFace.WallVelocity);
            if (wallFace.Contact) {
                wallSeparationKktResidual = std::max(
                    wallSeparationKktResidual,
                    std::max(-contactPressure(wallFace), 0.0f));
            } else {
                separatingCell[wallFace.Row] = 1;
                wallSeparationKktResidual = std::max(
                    wallSeparationKktResidual,
                    std::max(-separationSpeed, 0.0f));
                if (separationSpeed > 1e-5f)
                    ++separatingWallFaces;
            }
        }

        for (int row = 0; row < int(rowToCell.size()); ++row) {
            int const idx = rowToCell[row];
            m_p[idx]      = float(pressure[row]);
            separatingWallCells += separatingCell[row] ? 1 : 0;
        }
    }

    void FreeSurfaceSeparationSimulator::updateWallContactDiagnostics() {
        float const leftWall        = -0.5f + m_h;
        float const contactDistance = 2.0f * m_h;
        wallContactParticles        = 0;
        averageLeftWallDistance     = 0.0f;

        for (int particle = 0; particle < m_iNumSpheres; ++particle) {
            float const distance = std::max(m_particlePos[particle].x - leftWall, 0.0f);
            averageLeftWallDistance += distance;
            if (distance <= contactDistance) {
                ++wallContactParticles;
                if (highlightWallContact && m_colorMode == ColorMode::Default)
                    m_particleColor[particle] = glm::vec3(1.0f, 0.35f, 0.08f);
            }
        }

        if (m_iNumSpheres > 0)
            averageLeftWallDistance /= float(m_iNumSpheres);
    }
} // namespace VCX::Labs::Final
