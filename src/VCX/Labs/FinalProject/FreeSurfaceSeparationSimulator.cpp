#include "Labs/FinalProject/FreeSurfaceSeparationSimulator.h"

#include <array>
#include <cmath>

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
            FaceNeighbor { {  1, 0, 0 }, { 1, 0, 0 }, 0,  1.0f },
            FaceNeighbor { { 0, -1, 0 }, { 0, 0, 0 }, 1, -1.0f },
            FaceNeighbor { { 0,  1, 0 }, { 0, 1, 0 }, 1,  1.0f },
            FaceNeighbor { { 0, 0, -1 }, { 0, 0, 0 }, 2, -1.0f },
            FaceNeighbor { { 0, 0,  1 }, { 0, 0, 1 }, 2,  1.0f },
        };
    }

    void FreeSurfaceSeparationSimulator::setupScene(int res) {
        Simulator::setupScene(res);

        m_body = nullptr;
        gravity = { 0.0f, -9.81f, 0.0f };
        m_fRatio = 0.97f;
        numSubSteps = 2;
        numParticleIters = 3;
        numPressureIters = 120;
        separateParticles = true;
        compensateDrift = false;

        setupSplashParticles();
        wallSeparationCandidates = 0;
        wallSeparationClampedCells = 0;
        wallSeparationActiveSetIterations = 0;
        separatingWallFaces = 0;
        wallContactParticles = 0;
        minimumUnconstrainedPressure = 0.0f;
        pressureResidual = 0.0f;
        averageLeftWallDistance = 0.0f;
        updateWallContactDiagnostics();
    }

    void FreeSurfaceSeparationSimulator::setupSplashParticles() {
        // A thin sheet makes the wall attachment visible from the default
        // front camera, matching the 2D comparison in Figure 6 of the paper.
        glm::vec3 const center { 0.12f, 0.16f, 0.0f };
        glm::vec3 const radius { 0.18f, 0.18f, 0.08f };
        float const spacing = std::max(2.0f * m_particleRadius, 0.42f * m_h);

        std::vector<glm::vec3> positions;
        for (float x = center.x - radius.x; x <= center.x + radius.x; x += spacing) {
            int const xLayer = int(std::round((x - center.x + radius.x) / spacing));
            for (float y = center.y - radius.y; y <= center.y + radius.y; y += spacing) {
                int const yLayer = int(std::round((y - center.y + radius.y) / spacing));
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
        bool touchesSolid = false;
        bool touchesAir = false;

        for (FaceNeighbor const & side : FaceNeighbors) {
            glm::ivec3 const neighbor = cell + side.CellOffset;
            if (isSolidCell(neighbor)) {
                touchesSolid = true;
            } else if (m_type[gridOffset(neighbor)] == EMPTY_CELL) {
                touchesAir = true;
            }
        }
        return touchesSolid && touchesAir;
    }

    void FreeSurfaceSeparationSimulator::solveIncompressibility(
        int numIters,
        float dt,
        float overRelaxation,
        bool compensateDrift) {
        (void) dt;
        (void) overRelaxation;
        (void) compensateDrift;

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
        wallSeparationCandidates = 0;
        wallSeparationClampedCells = 0;
        wallSeparationActiveSetIterations = 0;
        separatingWallFaces = 0;
        minimumUnconstrainedPressure = 0.0f;
        pressureResidual = 0.0f;
        if (rowToCell.empty())
            return;

        auto decodeCell = [&](int idx) {
            return glm::ivec3(
                idx % m_iCellX,
                (idx / m_iCellX) % m_iCellY,
                idx / (m_iCellX * m_iCellY));
        };

        std::vector<char> candidate(rowToCell.size(), 0);
        std::vector<char> activeDirichlet(rowToCell.size(), 0);
        for (int row = 0; row < int(rowToCell.size()); ++row) {
            candidate[row] = isWallSeparationCandidate(decodeCell(rowToCell[row])) ? 1 : 0;
            wallSeparationCandidates += candidate[row] ? 1 : 0;
        }

        Eigen::VectorXd pressure = Eigen::VectorXd::Zero(int(rowToCell.size()));
        auto solveWithActiveSet = [&]() {
            int const matrixSize = int(rowToCell.size());
            Eigen::SparseMatrix<double> matrix(matrixSize, matrixSize);
            std::vector<Eigen::Triplet<double>> triplets;
            triplets.reserve(matrixSize * 7);
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(matrixSize);

            for (int row = 0; row < matrixSize; ++row) {
                if (activeDirichlet[row]) {
                    triplets.emplace_back(row, row, 1.0);
                    continue;
                }

                glm::ivec3 const cell = decodeCell(rowToCell[row]);
                double diagonal = 0.0;
                double divergence = 0.0;
                for (FaceNeighbor const & side : FaceNeighbors) {
                    glm::ivec3 const face = cell + side.FaceOffset;
                    glm::ivec3 const neighbor = cell + side.CellOffset;
                    int const faceIdx = gridOffset(face);
                    divergence += side.DivergenceSign * double(m_vel[faceIdx][side.Direction]);
                    diagonal += 1.0;

                    if (! isSolidCell(neighbor) && m_type[gridOffset(neighbor)] == FLUID_CELL) {
                        int const neighborRow = cellToRow[gridOffset(neighbor)];
                        if (neighborRow >= 0 && ! activeDirichlet[neighborRow])
                            triplets.emplace_back(row, neighborRow, -1.0);
                    }
                }

                triplets.emplace_back(row, row, diagonal);
                rhs[row] = -divergence;
            }

            matrix.setFromTriplets(triplets.begin(), triplets.end());
            Eigen::ConjugateGradient<
                Eigen::SparseMatrix<double>,
                Eigen::Lower | Eigen::Upper,
                Eigen::IncompleteCholesky<double>> solver;
            solver.setMaxIterations(std::max(numIters, 1));
            solver.setTolerance(1e-5);
            solver.compute(matrix);
            pressure = solver.solve(rhs);
            pressureResidual = float(solver.error());

            for (int row = 0; row < int(activeDirichlet.size()); ++row) {
                if (activeDirichlet[row])
                    pressure[row] = 0.0;
            }
        };

        bool changedActiveSet = false;
        int const maxActiveSetIterations = enableWallSeparation
            ? std::min<int>(32, int(rowToCell.size()) + 1)
            : 1;
        for (int iter = 0; iter < maxActiveSetIterations; ++iter) {
            solveWithActiveSet();
            wallSeparationActiveSetIterations = iter + 1;

            if (iter == 0 && pressure.size() > 0)
                minimumUnconstrainedPressure = float(pressure.minCoeff());

            changedActiveSet = false;
            if (! enableWallSeparation)
                break;

            for (int row = 0; row < int(rowToCell.size()); ++row) {
                if (candidate[row] && ! activeDirichlet[row] && pressure[row] < -1e-5) {
                    activeDirichlet[row] = 1;
                    changedActiveSet = true;
                }
            }
            if (! changedActiveSet)
                break;
        }

        if (changedActiveSet) {
            solveWithActiveSet();
            ++wallSeparationActiveSetIterations;
        }

        for (int row = 0; row < int(rowToCell.size()); ++row) {
            int const idx = rowToCell[row];
            m_p[idx] = float(pressure[row]);
            wallSeparationClampedCells += activeDirichlet[row] ? 1 : 0;
        }

        for (int row = 0; row < int(rowToCell.size()); ++row) {
            glm::ivec3 const cell = decodeCell(rowToCell[row]);
            float const cellPressure = float(pressure[row]);
            for (FaceNeighbor const & side : FaceNeighbors) {
                glm::ivec3 const face = cell + side.FaceOffset;
                int const faceIdx = gridOffset(face);
                m_vel[faceIdx][side.Direction] += side.DivergenceSign * cellPressure;
            }
        }

        enforceWallVelocityCondition();
    }

    void FreeSurfaceSeparationSimulator::enforceWallVelocityCondition() {
        separatingWallFaces = 0;
        for (int idx = 0; idx < m_iNumCells; ++idx) {
            if (m_type[idx] != FLUID_CELL)
                continue;

            glm::ivec3 const cell {
                idx % m_iCellX,
                (idx / m_iCellX) % m_iCellY,
                idx / (m_iCellX * m_iCellY),
            };

            for (FaceNeighbor const & side : FaceNeighbors) {
                glm::ivec3 const neighbor = cell + side.CellOffset;
                if (! isSolidCell(neighbor))
                    continue;

                glm::ivec3 const face = cell + side.FaceOffset;
                float & velocity = m_vel[gridOffset(face)][side.Direction];
                if (! enableWallSeparation) {
                    velocity = 0.0f;
                    continue;
                }

                // The solid normal points into the fluid. Negative relative normal
                // velocity penetrates the wall; positive velocity is free separation.
                if (side.DivergenceSign < 0.0f)
                    velocity = std::max(velocity, 0.0f);
                else
                    velocity = std::min(velocity, 0.0f);

                if (std::abs(velocity) > 1e-5f)
                    ++separatingWallFaces;
            }
        }
    }

    void FreeSurfaceSeparationSimulator::updateWallContactDiagnostics() {
        float const leftWall = -0.5f + m_h;
        float const contactDistance = 2.0f * m_h;
        wallContactParticles = 0;
        averageLeftWallDistance = 0.0f;

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
}
