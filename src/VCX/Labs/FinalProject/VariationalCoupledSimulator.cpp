#include "Labs/FinalProject/VariationalCoupledSimulator.h"

#include <algorithm>
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
    }

    void VariationalCoupledSimulator::setupScene(int res) {
        Simulator::setupScene(res);
        numPressureIters = 120;
        for (auto & fractions : _faceFluidFraction)
            fractions.assign(m_iNumCells, 0.0f);
        pressureResidual       = 0.0f;
        pressureSolveSucceeded = true;
    }

    bool VariationalCoupledSimulator::isValidCell(glm::ivec3 const & cell) const {
        return cell.x >= 0 && cell.x < m_iCellX
            && cell.y >= 0 && cell.y < m_iCellY
            && cell.z >= 0 && cell.z < m_iCellZ;
    }

    int VariationalCoupledSimulator::gridOffset(glm::ivec3 const & cell) const {
        return cell.x + m_iCellX * (cell.y + m_iCellY * cell.z);
    }

    glm::vec3 VariationalCoupledSimulator::pressureCellCenter(glm::ivec3 const & cell) const {
        return (glm::vec3(cell) + glm::vec3(0.5f)) * m_h - glm::vec3(0.5f);
    }

    glm::vec3 VariationalCoupledSimulator::faceCenter(glm::ivec3 const & face, int dir) const {
        glm::vec3 offset(0.5f);
        offset[dir] = 0.0f;
        return (glm::vec3(face) + offset) * m_h - glm::vec3(0.5f);
    }

    bool VariationalCoupledSimulator::isTankFaceOpen(glm::ivec3 const & face, int dir) const {
        glm::ivec3 left = face;
        left[dir] -= 1;
        glm::ivec3 const right = face;
        if (! isValidCell(left) || ! isValidCell(right))
            return false;
        return m_s[gridOffset(left)] > 0.0f
            && m_s[gridOffset(right)] > 0.0f;
    }

    bool VariationalCoupledSimulator::isSolidPressureCell(glm::ivec3 const & cell) const {
        if (! isValidCell(cell) || m_s[gridOffset(cell)] <= 0.0f)
            return true;
        return m_body && m_body->GetSDF(pressureCellCenter(cell)) < -1e-5f;
    }

    bool VariationalCoupledSimulator::isWallSeparationCandidate(glm::ivec3 const & cell) const {
        bool touchesSolid = false;
        bool touchesAir   = false;

        for (FaceNeighbor const & side : FaceNeighbors) {
            glm::ivec3 const neighbor = cell + side.CellOffset;
            if (isSolidPressureCell(neighbor)) {
                touchesSolid = true;
            } else if (m_type[gridOffset(neighbor)] == EMPTY_CELL) {
                touchesAir = true;
            }
        }
        return touchesSolid && touchesAir;
    }

    float VariationalCoupledSimulator::estimateFaceFluidFraction(
        glm::ivec3 const & face,
        int                dir) const {
        if (! isTankFaceOpen(face, dir))
            return 0.0f;
        if (! m_body)
            return 1.0f;

        glm::vec3 const center         = faceCenter(face, dir);
        float const     centerDistance = m_body->GetSDF(center);
        float const     sampleRadius   = 0.5f * std::sqrt(3.0f) * m_h;
        if (centerDistance >= sampleRadius)
            return 1.0f;
        if (centerDistance <= -sampleRadius)
            return 0.0f;

        int const sampleCount  = std::clamp(volumeSamplesPerAxis, 1, 16);
        int       fluidSamples = 0;

        for (int x = 0; x < sampleCount; ++x) {
            for (int y = 0; y < sampleCount; ++y) {
                for (int z = 0; z < sampleCount; ++z) {
                    glm::vec3 const unitSample {
                        (float(x) + 0.5f) / float(sampleCount) - 0.5f,
                        (float(y) + 0.5f) / float(sampleCount) - 0.5f,
                        (float(z) + 0.5f) / float(sampleCount) - 0.5f,
                    };
                    if (m_body->GetSDF(center + unitSample * m_h) >= 0.0f)
                        ++fluidSamples;
                }
            }
        }

        return float(fluidSamples)
            / float(sampleCount * sampleCount * sampleCount);
    }

    void VariationalCoupledSimulator::updateFaceFluidFractions() {
        for (int k = 0; k < m_iCellZ; ++k) {
            for (int j = 0; j < m_iCellY; ++j) {
                for (int i = 0; i < m_iCellX; ++i) {
                    glm::ivec3 const face { i, j, k };
                    int const        idx = gridOffset(face);
                    for (int dir = 0; dir < 3; ++dir)
                        _faceFluidFraction[dir][idx] = estimateFaceFluidFraction(face, dir);
                }
            }
        }
    }

    float VariationalCoupledSimulator::faceWeight(int dir, int faceIdx) const {
        float const fraction = _faceFluidFraction[dir][faceIdx];
        if (useSubgridWeights)
            return fraction;
        return fraction > 0.5f ? 1.0f : 0.0f;
    }

    float VariationalCoupledSimulator::solidVelocity(
        glm::ivec3 const & face,
        int                dir,
        glm::ivec3 const & solidCell) const {
        if (! m_body || ! isValidCell(solidCell))
            return 0.0f;

        glm::vec3 const center = faceCenter(face, dir);
        if (m_body->GetSDF(pressureCellCenter(solidCell)) >= 0.0f
            && m_body->GetSDF(center) >= 0.0f)
            return 0.0f;

        return m_body->GetVelocityAtPoint(center - m_body->position)[dir];
    }

    void VariationalCoupledSimulator::solveIncompressibility(
        int   numIters,
        float dt,
        float overRelaxation,
        bool  compensateDrift) {
        (void) overRelaxation;

        updateFaceFluidFractions();

        std::vector<int> cellToRow(m_iNumCells, -1);
        std::vector<int> rowToCell;
        rowToCell.reserve(m_iNumCells / 3);
        for (int idx = 0; idx < m_iNumCells; ++idx) {
            glm::ivec3 const cell {
                idx % m_iCellX,
                (idx / m_iCellX) % m_iCellY,
                idx / (m_iCellX * m_iCellY),
            };
            if (m_type[idx] == FLUID_CELL && ! isSolidPressureCell(cell)) {
                cellToRow[idx] = int(rowToCell.size());
                rowToCell.push_back(idx);
            }
        }

        std::fill(m_p.begin(), m_p.end(), 0.0f);
        m_feedbackForce        = glm::vec3(0.0f);
        m_feedbackTorque       = glm::vec3(0.0f);
        pressureResidual       = 0.0f;
        pressureSolveSucceeded = true;
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
            float BoundaryWeight;
            float WallVelocity;
            bool  Candidate;
            bool  Contact;
        };

        std::vector<WallFace> wallFaces;
        std::vector<int>      wallFaceForSide(rowToCell.size() * FaceNeighbors.size(), -1);
        for (int row = 0; row < int(rowToCell.size()); ++row) {
            glm::ivec3 const cell = decodeCell(rowToCell[row]);
            bool const       candidateCell =
                enableWallSeparation && isWallSeparationCandidate(cell);
            for (int sideIndex = 0; sideIndex < int(FaceNeighbors.size()); ++sideIndex) {
                FaceNeighbor const & side     = FaceNeighbors[sideIndex];
                glm::ivec3 const     neighbor = cell + side.CellOffset;
                if (! isSolidPressureCell(neighbor))
                    continue;

                glm::ivec3 const face    = cell + side.FaceOffset;
                int const        faceIdx = gridOffset(face);
                float const      boundaryWeight =
                    1.0f - faceWeight(side.Direction, faceIdx);
                if (boundaryWeight <= 1e-6f)
                    continue;

                wallFaceForSide[row * FaceNeighbors.size() + sideIndex] =
                    int(wallFaces.size());
                wallFaces.push_back(WallFace {
                    .Row            = row,
                    .FaceIndex      = faceIdx,
                    .Direction      = side.Direction,
                    .DivergenceSign = side.DivergenceSign,
                    .BoundaryWeight = boundaryWeight,
                    .WallVelocity   = solidVelocity(face, side.Direction, neighbor),
                    .Candidate      = candidateCell,
                    .Contact        = true,
                });
            }
        }

        std::vector<glm::vec3> const intermediateVelocity = m_vel;
        Eigen::VectorXd              pressure             = Eigen::VectorXd::Zero(int(rowToCell.size()));
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
                        int const            faceIdx  = gridOffset(cell + side.FaceOffset);
                        int const            wallFaceIndex =
                            wallFaceForSide[row * FaceNeighbors.size() + sideIndex];
                        if (wallFaceIndex >= 0) {
                            WallFace const & wallFace = wallFaces[wallFaceIndex];
                            if (wallFace.Candidate && ! wallFace.Contact)
                                hasDirichletBoundary = true;
                            continue;
                        }
                        if (isSolidPressureCell(neighbor))
                            continue;

                        float const openWeight = faceWeight(side.Direction, faceIdx);
                        if (openWeight <= 1e-6f)
                            continue;

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

                int const        idx                = rowToCell[row];
                glm::ivec3 const cell               = decodeCell(idx);
                double           diagonal           = 0.0;
                double           weightedDivergence = 0.0;

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
                        weightedDivergence += side.DivergenceSign
                            * double(wallFace.BoundaryWeight) * double(velocity);
                        if (wallFace.Candidate && ! wallFace.Contact)
                            diagonal += wallFace.BoundaryWeight;
                        continue;
                    }

                    if (isSolidPressureCell(neighbor))
                        continue;

                    float const openWeight = faceWeight(side.Direction, faceIdx);
                    if (openWeight <= 1e-6f)
                        continue;

                    weightedDivergence += side.DivergenceSign
                        * double(openWeight)
                        * double(intermediateVelocity[faceIdx][side.Direction]);
                    diagonal += openWeight;

                    int const neighborRow = cellToRow[gridOffset(neighbor)];
                    if (neighborRow >= 0 && ! pinnedRows[neighborRow])
                        triplets.emplace_back(row, neighborRow, -double(openWeight));
                }

                if (compensateDrift && m_particleRestDensity > 0.0f) {
                    float const densityError = m_particleDensity[idx] - m_particleRestDensity;
                    if (densityError > 0.0f)
                        weightedDivergence -= compensateDriftWeight * densityError;
                }

                if (diagonal <= 1e-8) {
                    diagonal           = 1.0;
                    weightedDivergence = 0.0;
                }
                triplets.emplace_back(row, row, diagonal);
                rhs[row] = -weightedDivergence;
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
            pressureSolveSucceeded =
                solver.info() == Eigen::Success && pressure.allFinite();
            if (pressureSolveSucceeded) {
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
                glm::ivec3 const cell          = decodeCell(rowToCell[row]);
                float const      pressureValue = float(pressure[row]);
                for (int sideIndex = 0; sideIndex < int(FaceNeighbors.size()); ++sideIndex) {
                    FaceNeighbor const & side     = FaceNeighbors[sideIndex];
                    glm::ivec3 const     face     = cell + side.FaceOffset;
                    glm::ivec3 const     neighbor = cell + side.CellOffset;
                    int const            faceIdx  = gridOffset(face);
                    int const            wallFaceIndex =
                        wallFaceForSide[row * FaceNeighbors.size() + sideIndex];
                    if (wallFaceIndex >= 0) {
                        WallFace const & wallFace = wallFaces[wallFaceIndex];
                        if (! wallFace.Candidate || wallFace.Contact)
                            continue;
                    } else if (
                        isSolidPressureCell(neighbor)
                        || faceWeight(side.Direction, faceIdx) <= 1e-6f) {
                        continue;
                    }

                    m_vel[faceIdx][side.Direction] +=
                        side.DivergenceSign * pressureValue;
                }
            }
        };

        bool      changedActiveSet       = false;
        int const maxActiveSetIterations = enableWallSeparation
            ? std::min<int>(256, 2 * int(wallFaces.size()) + 8)
            : 1;
        for (int iter = 0; iter < maxActiveSetIterations; ++iter) {
            solveWithContacts();
            changedActiveSet = false;
            if (! enableWallSeparation)
                break;

            bool removedNegativePressureContact = false;
            for (WallFace & wallFace : wallFaces) {
                if (! wallFace.Candidate || ! wallFace.Contact)
                    continue;

                if (pressure[wallFace.Row] < -1e-5) {
                    wallFace.Contact = false;
                    removedNegativePressureContact = true;
                    changedActiveSet = true;
                }
            }
            if (removedNegativePressureContact)
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
        if (changedActiveSet)
            solveWithContacts();

        for (int row = 0; row < int(rowToCell.size()); ++row)
            m_p[rowToCell[row]] = float(pressure[row]);

        std::vector<BoundaryPressureSample> boundarySamples;
        boundarySamples.reserve(wallFaces.size());
        for (WallFace const & wallFace : wallFaces) {
            boundarySamples.push_back(BoundaryPressureSample {
                .PressureCell    = rowToCell[wallFace.Row],
                .FaceIndex       = wallFace.FaceIndex,
                .Direction       = wallFace.Direction,
                .DivergenceSign  = wallFace.DivergenceSign,
                .BoundaryWeight  = wallFace.BoundaryWeight,
                .Contact         = ! wallFace.Candidate || wallFace.Contact,
            });
        }

        applyRigidBodyFeedback(boundarySamples, dt);
    }

    void VariationalCoupledSimulator::applyRigidBodyFeedback(
        std::vector<BoundaryPressureSample> const & boundarySamples,
        float                                      dt) {
        if (! m_body)
            return;

        float const invDt = 1.0f / std::max(dt, 1e-6f);
        for (BoundaryPressureSample const & sample : boundarySamples) {
            if (! sample.Contact || sample.BoundaryWeight <= 1e-6f)
                continue;

            float const pressure = m_p[sample.PressureCell];
            if (pressure <= 0.0f)
                continue;

            glm::ivec3 const face {
                sample.FaceIndex % m_iCellX,
                (sample.FaceIndex / m_iCellX) % m_iCellY,
                sample.FaceIndex / (m_iCellX * m_iCellY),
            };
            glm::vec3 const applicationPoint = faceCenter(face, sample.Direction);
            glm::vec3       forceDirection(0.0f);
            forceDirection[sample.Direction] = sample.DivergenceSign;
            float const     area  = sample.BoundaryWeight * m_h * m_h;
            glm::vec3 const force = forceDirection * pressure * area * invDt;
            m_feedbackForce += force;
            m_feedbackTorque += glm::cross(applicationPoint - m_body->position, force);
        }
    }
} // namespace VCX::Labs::Final
