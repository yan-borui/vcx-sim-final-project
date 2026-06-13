#include "Labs/FinalProject/VariationalCoupledSimulator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

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
            FaceNeighbor { { 0, 0, 1 }, { 0, 0, 1 }, 2,  1.0f },
        };

        using Matrix6 = Eigen::Matrix<double, 6, 6>;

        struct BoundQpResult {
            Eigen::VectorXd Pressure;
            Eigen::VectorXd Gradient;
            int             Sweeps;
            bool            Converged;
        };

        BoundQpResult solveNonnegativeQuadratic(
            Eigen::MatrixXd const & hessian,
            Eigen::VectorXd const & rhs) {
            int const size = int(rhs.size());
            BoundQpResult result {
                .Pressure  = Eigen::VectorXd::Zero(size),
                .Gradient  = -rhs,
                .Sweeps    = 0,
                .Converged = size == 0,
            };
            if (size == 0)
                return result;

            constexpr double KktTolerance = 1e-7;
            std::vector<char> freeVariables(size, 0);
            for (int index = 0; index < size; ++index)
                freeVariables[index] = rhs[index] > KktTolerance;

            // Solve the current contact set exactly. A repeated set indicates
            // degeneracy, in which case projected coordinate descent finishes
            // the same convex QP and is accepted only after a full KKT check.
            std::unordered_set<std::size_t> visitedSets;
            for (int iteration = 0; iteration < 4 * size + 32; ++iteration) {
                std::size_t signature = 1469598103934665603ull;
                for (char const isFree : freeVariables) {
                    signature ^= std::size_t(isFree);
                    signature *= 1099511628211ull;
                }
                if (! visitedSets.insert(signature).second)
                    break;

                std::vector<int> freeIndices;
                freeIndices.reserve(size);
                for (int index = 0; index < size; ++index) {
                    if (freeVariables[index])
                        freeIndices.push_back(index);
                }

                result.Pressure.setZero();
                if (! freeIndices.empty()) {
                    int const freeCount = int(freeIndices.size());
                    Eigen::MatrixXd freeHessian(freeCount, freeCount);
                    Eigen::VectorXd freeRhs(freeCount);
                    for (int row = 0; row < freeCount; ++row) {
                        freeRhs[row] = rhs[freeIndices[row]];
                        for (int column = 0; column < freeCount; ++column) {
                            freeHessian(row, column) =
                                hessian(
                                    freeIndices[row],
                                    freeIndices[column]);
                        }
                    }
                    Eigen::LDLT<Eigen::MatrixXd> freeFactor(freeHessian);
                    if (freeFactor.info() != Eigen::Success)
                        break;
                    Eigen::VectorXd const freePressure =
                        freeFactor.solve(freeRhs);
                    if (! freePressure.allFinite())
                        break;
                    for (int row = 0; row < freeCount; ++row)
                        result.Pressure[freeIndices[row]] =
                            freePressure[row];
                }

                bool removedNegativePressure = false;
                for (int index = 0; index < size; ++index) {
                    if (freeVariables[index]
                        && result.Pressure[index] < -KktTolerance) {
                        freeVariables[index] = 0;
                        removedNegativePressure = true;
                    }
                }
                if (removedNegativePressure)
                    continue;

                result.Pressure = result.Pressure.cwiseMax(0.0);
                result.Gradient =
                    hessian * result.Pressure - rhs;

                bool releasedConstraint = false;
                for (int index = 0; index < size; ++index) {
                    if (! freeVariables[index]
                        && result.Gradient[index] < -KktTolerance) {
                        freeVariables[index] = 1;
                        releasedConstraint = true;
                    }
                }
                result.Sweeps = iteration + 1;
                if (! releasedConstraint) {
                    result.Converged = true;
                    return result;
                }
            }

            result.Pressure = result.Pressure.cwiseMax(0.0);
            result.Gradient = hessian * result.Pressure - rhs;
            constexpr int MaximumSweeps = 100000;
            for (int sweep = 0; sweep < MaximumSweeps; ++sweep) {
                for (int index = 0; index < size; ++index) {
                    double const diagonal = hessian(index, index);
                    if (diagonal <= 1e-12)
                        continue;
                    double const nextPressure = std::max(
                        0.0,
                        result.Pressure[index]
                            - result.Gradient[index] / diagonal);
                    double const delta =
                        nextPressure - result.Pressure[index];
                    if (std::abs(delta) <= 1e-15)
                        continue;
                    result.Pressure[index] = nextPressure;
                    result.Gradient.noalias() +=
                        delta * hessian.col(index);
                }

                double kktResidual = 0.0;
                for (int index = 0; index < size; ++index) {
                    double const violation =
                        result.Pressure[index] > KktTolerance
                        ? std::abs(result.Gradient[index])
                        : std::max(-result.Gradient[index], 0.0);
                    kktResidual = std::max(kktResidual, violation);
                }
                ++result.Sweeps;
                if (kktResidual <= KktTolerance) {
                    result.Converged = true;
                    break;
                }
            }
            return result;
        }
    }

    void VariationalCoupledSimulator::setupScene(int res) {
        Simulator::setupScene(res);
        voxelizeDynamicBody = false;
        numPressureIters = 120;
        for (auto & fractions : _faceOpenFraction)
            fractions.assign(m_iNumCells, 0.0f);
        for (auto & fractions : _faceFluidFraction)
            fractions.assign(m_iNumCells, 0.0f);
        pressureResidual            = 0.0f;
        wallSeparationKktResidual   = 0.0f;
        wallSeparationQpResidual    = 0.0f;
        wallSeparationSpeedResidual = 0.0f;
        maximumBoundaryPressure     = 0.0f;
        wallSeparationIterations    = 0;
        wallSeparationCandidateCount = 0;
        _wallPressureWarmStart.clear();
        pressureSolveSucceeded      = true;
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

    bool VariationalCoupledSimulator::hasFluidSupportAcrossOpenFace(
        glm::ivec3 const & cell) const {
        if (! isSolidPressureCell(cell)
            || ! isValidCell(cell)
            || m_s[gridOffset(cell)] <= 0.0f)
            return false;

        for (FaceNeighbor const & side : FaceNeighbors) {
            glm::ivec3 const face = cell + side.FaceOffset;
            if (! isValidCell(face))
                continue;
            int const faceIdx = gridOffset(face);
            if (faceWeight(side.Direction, faceIdx) <= 1e-6f)
                continue;

            glm::ivec3 const neighbor = cell + side.CellOffset;
            if (! isValidCell(neighbor) || m_s[gridOffset(neighbor)] <= 0.0f)
                continue;
            if (m_type[gridOffset(neighbor)] == FLUID_CELL)
                return true;
        }
        return false;
    }

    bool VariationalCoupledSimulator::isPressureUnknownCell(
        glm::ivec3 const & cell) const {
        if (! isValidCell(cell) || m_s[gridOffset(cell)] <= 0.0f)
            return false;
        return m_type[gridOffset(cell)] == FLUID_CELL
            || hasFluidSupportAcrossOpenFace(cell);
    }

    bool VariationalCoupledSimulator::isWallSeparationCandidate(glm::ivec3 const & cell) const {
        for (FaceNeighbor const & side : FaceNeighbors) {
            glm::ivec3 const face     = cell + side.FaceOffset;
            int const        faceIdx  = gridOffset(face);
            bool const       tankOpen = isTankFaceOpen(face, side.Direction);
            float const      boundaryWeight =
                tankOpen ? 1.0f - _faceOpenFraction[side.Direction][faceIdx] : 1.0f;
            if (boundaryWeight > 1e-6f)
                return true;
        }
        return false;
    }

    glm::vec2 VariationalCoupledSimulator::estimateFaceFractions(
        glm::ivec3 const &         face,
        int                        dir,
        std::vector<float> const & liquidLevelSet) const {
        if (! isTankFaceOpen(face, dir))
            return glm::vec2(0.0f);

        glm::vec3 const center      = faceCenter(face, dir);
        int const sampleCount       = std::clamp(volumeSamplesPerAxis, 1, 16);
        int       openSamples       = 0;
        int       fluidSamples      = 0;
        for (int x = 0; x < sampleCount; ++x) {
            for (int y = 0; y < sampleCount; ++y) {
                for (int z = 0; z < sampleCount; ++z) {
                    glm::vec3 const unitSample {
                        (float(x) + 0.5f) / float(sampleCount) - 0.5f,
                        (float(y) + 0.5f) / float(sampleCount) - 0.5f,
                        (float(z) + 0.5f) / float(sampleCount) - 0.5f,
                    };
                    glm::vec3 const samplePosition =
                        center + unitSample * m_h;
                    if (m_body && m_body->GetSDF(samplePosition) < 0.0f)
                        continue;

                    ++openSamples;
                    if (sampleCellCenteredField(
                            liquidLevelSet,
                            samplePosition)
                        <= 0.0f)
                        ++fluidSamples;
                }
            }
        }
        float const inverseSampleCount =
            1.0f / float(sampleCount * sampleCount * sampleCount);
        return {
            float(openSamples) * inverseSampleCount,
            float(fluidSamples) * inverseSampleCount,
        };
    }

    void VariationalCoupledSimulator::updateFaceFluidFractions(
        std::vector<float> const & liquidLevelSet) {
        for (int k = 0; k < m_iCellZ; ++k) {
            for (int j = 0; j < m_iCellY; ++j) {
                for (int i = 0; i < m_iCellX; ++i) {
                    glm::ivec3 const face { i, j, k };
                    int const        idx = gridOffset(face);
                    for (int dir = 0; dir < 3; ++dir) {
                        glm::vec2 const fractions =
                            estimateFaceFractions(face, dir, liquidLevelSet);
                        _faceOpenFraction[dir][idx]  = fractions.x;
                        _faceFluidFraction[dir][idx] = fractions.y;
                    }
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

    float VariationalCoupledSimulator::wallGhostPressureScale(
        glm::ivec3 const & cell,
        glm::ivec3 const & face,
        int                dir,
        bool               bodyBoundary) const {
        if (! bodyBoundary || ! m_body)
            return 2.0f;

        float const minimumDistance = 1e-4f * m_h;
        glm::vec3 const cellCenter  = pressureCellCenter(cell);
        glm::vec3 const wallSample  = faceCenter(face, dir);
        float const     cellPhi =
            std::max(m_body->GetSDF(cellCenter), minimumDistance);
        float const samplePhi = m_body->GetSDF(wallSample);

        float distance = cellPhi;
        if (samplePhi < -minimumDistance) {
            float const segmentFraction = std::clamp(
                cellPhi / (cellPhi - samplePhi),
                0.05f,
                1.0f);
            distance =
                segmentFraction * glm::length(wallSample - cellCenter);
        }
        return m_h / std::clamp(distance, 0.05f * m_h, m_h);
    }

    float VariationalCoupledSimulator::solidVelocity(
        glm::ivec3 const & face,
        int                dir,
        glm::ivec3 const & solidCell) const {
        if (! m_body || ! isValidCell(solidCell))
            return 0.0f;
        glm::vec3 const center = faceCenter(face, dir);
        return m_body->GetVelocityAtPoint(center - m_body->position)[dir];
    }

    void VariationalCoupledSimulator::solveIncompressibility(
        int   numIters,
        float dt,
        float overRelaxation,
        bool  compensateDrift) {
        (void) numIters;
        (void) overRelaxation;

        std::vector<float> const liquidLevelSet = buildParticleLevelSet();
        updateFaceFluidFractions(liquidLevelSet);

        std::vector<int> cellToRow(m_iNumCells, -1);
        std::vector<int> rowToCell;
        rowToCell.reserve(m_iNumCells / 3);
        for (int idx = 0; idx < m_iNumCells; ++idx) {
            glm::ivec3 const cell {
                idx % m_iCellX,
                (idx / m_iCellX) % m_iCellY,
                idx / (m_iCellX * m_iCellY),
            };
            if (isPressureUnknownCell(cell)) {
                cellToRow[idx] = int(rowToCell.size());
                rowToCell.push_back(idx);
            }
        }

        std::fill(m_p.begin(), m_p.end(), 0.0f);
        m_feedbackForce             = glm::vec3(0.0f);
        m_feedbackTorque            = glm::vec3(0.0f);
        pressureResidual            = 0.0f;
        wallSeparationKktResidual   = 0.0f;
        wallSeparationIterations    = 0;
        pressureSolveSucceeded      = true;
        if (rowToCell.empty())
            return;

        auto decodeCell = [&](int idx) {
            return glm::ivec3(
                idx % m_iCellX,
                (idx / m_iCellX) % m_iCellY,
                idx / (m_iCellX * m_iCellY));
        };
        std::vector<glm::vec3> const intermediateVelocity = m_vel;

        struct WallFace {
            int   Row;
            int   FaceIndex;
            int   Direction;
            int   UnconstrainedIndex;
            int   CandidateIndex;
            float DivergenceSign;
            float OpenWeight;
            float BoundaryWeight;
            float GhostScale;
            float WallVelocity;
            bool  BodyBoundary;
            bool  Candidate;
        };

        int const pressureCount = int(rowToCell.size());
        std::vector<WallFace> wallFaces;
        for (int row = 0; row < pressureCount; ++row) {
            glm::ivec3 const cell = decodeCell(rowToCell[row]);
            bool const candidateCell =
                enableWallSeparation
                && m_type[rowToCell[row]] == FLUID_CELL
                && isWallSeparationCandidate(cell);
            for (int sideIndex = 0;
                 sideIndex < int(FaceNeighbors.size());
                 ++sideIndex) {
                FaceNeighbor const & side = FaceNeighbors[sideIndex];
                glm::ivec3 const neighbor = cell + side.CellOffset;
                glm::ivec3 const face     = cell + side.FaceOffset;
                int const        faceIdx  = gridOffset(face);
                bool const       tankOpen =
                    isTankFaceOpen(face, side.Direction);
                float const openWeight =
                    tankOpen ? faceWeight(side.Direction, faceIdx) : 0.0f;
                float const boundaryWeight =
                    tankOpen
                    ? 1.0f - _faceOpenFraction[side.Direction][faceIdx]
                    : 1.0f;
                if (boundaryWeight <= 1e-6f)
                    continue;

                bool const bodyBoundary = tankOpen && m_body;
                wallFaces.push_back(WallFace {
                    .Row                = row,
                    .FaceIndex          = faceIdx,
                    .Direction          = side.Direction,
                    .UnconstrainedIndex = -1,
                    .CandidateIndex     = -1,
                    .DivergenceSign     = side.DivergenceSign,
                    .OpenWeight         = openWeight,
                    .BoundaryWeight     = boundaryWeight,
                    .GhostScale         = wallGhostPressureScale(
                        cell,
                        face,
                        side.Direction,
                        bodyBoundary),
                    .WallVelocity       = solidVelocity(
                        face,
                        side.Direction,
                        neighbor),
                    .BodyBoundary       = bodyBoundary,
                    .Candidate          = candidateCell,
                });
            }
        }

        int unconstrainedCount = pressureCount;
        int candidateCount     = 0;
        for (WallFace & wallFace : wallFaces) {
            if (wallFace.Candidate)
                wallFace.CandidateIndex = candidateCount++;
            else
                wallFace.UnconstrainedIndex = unconstrainedCount++;
        }
        wallSeparationCandidateCount = candidateCount;
        std::vector<std::uint64_t> candidateKeys(candidateCount);
        for (WallFace const & wallFace : wallFaces) {
            if (! wallFace.Candidate)
                continue;
            candidateKeys[wallFace.CandidateIndex] =
                (std::uint64_t(std::uint32_t(wallFace.FaceIndex)) << 2)
                | std::uint64_t(wallFace.Direction);
        }

        std::vector<char> pinnedRows(pressureCount, 0);
        std::vector<char> visitedRows(pressureCount, 0);
        for (int seed = 0; seed < pressureCount; ++seed) {
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
                for (int sideIndex = 0;
                     sideIndex < int(FaceNeighbors.size());
                     ++sideIndex) {
                    FaceNeighbor const & side = FaceNeighbors[sideIndex];
                    glm::ivec3 const face = cell + side.FaceOffset;
                    glm::ivec3 const neighbor = cell + side.CellOffset;
                    int const faceIdx = gridOffset(face);
                    float const openWeight =
                        faceWeight(side.Direction, faceIdx);
                    if (openWeight <= 1e-6f
                        || ! isValidCell(neighbor)
                        || m_s[gridOffset(neighbor)] <= 0.0f)
                        continue;

                    int const neighborRow =
                        cellToRow[gridOffset(neighbor)];
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

        Eigen::VectorXd rhsY = Eigen::VectorXd::Zero(unconstrainedCount);
        Eigen::VectorXd rhsZ = Eigen::VectorXd::Zero(candidateCount);
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(unconstrainedCount * 9);
        std::vector<int> candidatePressureRows(candidateCount, -1);
        Eigen::VectorXd candidateCoupling =
            Eigen::VectorXd::Zero(candidateCount);
        Eigen::VectorXd candidateDiagonal =
            Eigen::VectorXd::Zero(candidateCount);
        Eigen::MatrixXd rigidY =
            Eigen::MatrixXd::Zero(unconstrainedCount, 6);
        Eigen::MatrixXd rigidZ =
            Eigen::MatrixXd::Zero(candidateCount, 6);

        for (WallFace const & wallFace : wallFaces) {
            int const   pressureRow = wallFace.Row;
            double const coefficient =
                double(wallFace.BoundaryWeight * wallFace.GhostScale);
            double const pressureRhs =
                -double(wallFace.DivergenceSign)
                * double(wallFace.BoundaryWeight)
                * double(intermediateVelocity[wallFace.FaceIndex][wallFace.Direction]);
            double const boundaryRhs =
                double(wallFace.DivergenceSign)
                * double(wallFace.BoundaryWeight)
                * double(
                    intermediateVelocity[wallFace.FaceIndex][wallFace.Direction]
                    - wallFace.WallVelocity);

            // q is pressure on the solid boundary. The fluid energy term is
            // c * (p - q)^2 / 2; only candidate q variables receive q >= 0.
            if (! pinnedRows[pressureRow]) {
                triplets.emplace_back(
                    pressureRow,
                    pressureRow,
                    coefficient);
                rhsY[pressureRow] += pressureRhs;
            }

            int qRow = wallFace.UnconstrainedIndex;
            if (wallFace.Candidate) {
                int const qColumn = wallFace.CandidateIndex;
                candidateDiagonal[qColumn] += coefficient;
                rhsZ[qColumn] += boundaryRhs;
                if (! pinnedRows[pressureRow]) {
                    candidatePressureRows[qColumn] = pressureRow;
                    candidateCoupling[qColumn] -= coefficient;
                }
            } else {
                triplets.emplace_back(qRow, qRow, coefficient);
                rhsY[qRow] += boundaryRhs;
                if (! pinnedRows[pressureRow]) {
                    triplets.emplace_back(
                        pressureRow,
                        qRow,
                        -coefficient);
                    triplets.emplace_back(
                        qRow,
                        pressureRow,
                        -coefficient);
                }
            }

            if (! wallFace.BodyBoundary
                || ! m_body
                || m_body->isStatic
                || m_body->mass <= 1e-6f)
                continue;

            glm::ivec3 const face = decodeCell(wallFace.FaceIndex);
            glm::vec3 direction(0.0f);
            direction[wallFace.Direction] = wallFace.DivergenceSign;
            glm::vec3 const linear =
                wallFace.BoundaryWeight * direction;
            glm::vec3 const lever =
                faceCenter(face, wallFace.Direction) - m_body->position;
            glm::vec3 const angular = glm::cross(lever, linear);

            auto assignRigidRow = [&](Eigen::MatrixXd & matrix, int row) {
                matrix(row, 0) = linear.x;
                matrix(row, 1) = linear.y;
                matrix(row, 2) = linear.z;
                matrix(row, 3) = angular.x;
                matrix(row, 4) = angular.y;
                matrix(row, 5) = angular.z;
            };
            if (wallFace.Candidate)
                assignRigidRow(rigidZ, wallFace.CandidateIndex);
            else
                assignRigidRow(rigidY, qRow);
        }

        for (int row = 0; row < pressureCount; ++row) {
            if (pinnedRows[row]) {
                triplets.emplace_back(row, row, 1.0);
                rhsY[row] = 0.0;
                continue;
            }

            int const        idx  = rowToCell[row];
            glm::ivec3 const cell = decodeCell(idx);
            double diagonal = 0.0;
            for (FaceNeighbor const & side : FaceNeighbors) {
                glm::ivec3 const face     = cell + side.FaceOffset;
                glm::ivec3 const neighbor = cell + side.CellOffset;
                int const        faceIdx  = gridOffset(face);
                float const openWeight =
                    faceWeight(side.Direction, faceIdx);
                if (openWeight <= 1e-6f
                    || ! isValidCell(neighbor)
                    || m_s[gridOffset(neighbor)] <= 0.0f)
                    continue;

                rhsY[row] -=
                    double(side.DivergenceSign)
                    * double(openWeight)
                    * double(intermediateVelocity[faceIdx][side.Direction]);

                int const neighborIdx = gridOffset(neighbor);
                int const neighborRow = cellToRow[neighborIdx];
                if (neighborRow >= 0) {
                    diagonal += openWeight;
                    if (! pinnedRows[neighborRow])
                        triplets.emplace_back(
                            row,
                            neighborRow,
                            -double(openWeight));
                } else if (m_type[neighborIdx] == EMPTY_CELL) {
                    diagonal += openWeight * ghostFluidPressureScale(
                        liquidLevelSet,
                        idx,
                        neighborIdx);
                } else {
                    diagonal += openWeight;
                }
            }

            if (compensateDrift && m_particleRestDensity > 0.0f) {
                float const densityError =
                    m_particleDensity[idx] - m_particleRestDensity;
                if (densityError > 0.0f)
                    rhsY[row] += compensateDriftWeight * densityError;
            }
            triplets.emplace_back(row, row, std::max(diagonal, 1e-8));
        }

        Eigen::SparseMatrix<double> baseYY(
            unconstrainedCount,
            unconstrainedCount);
        baseYY.setFromTriplets(triplets.begin(), triplets.end());

        bool const dynamicBody =
            m_body && ! m_body->isStatic && m_body->mass > 1e-6f;
        Eigen::MatrixXd lowRankY =
            Eigen::MatrixXd::Zero(unconstrainedCount, 6);
        Eigen::MatrixXd lowRankZ =
            Eigen::MatrixXd::Zero(candidateCount, 6);
        if (dynamicBody) {
            // Equation (13) adds J^T M_S^-1 J. Keep its rank-six factor
            // instead of inserting a dense pressure block.
            Matrix6 mobility = Matrix6::Zero();
            float const cellVolume = m_h * m_h * m_h;
            for (int axis = 0; axis < 3; ++axis)
                mobility(axis, axis) = cellVolume / m_body->mass;
            glm::mat3 const inertiaInv = m_body->GetInertiaWorldInv();
            for (int row = 0; row < 3; ++row) {
                for (int column = 0; column < 3; ++column) {
                    mobility(3 + row, 3 + column) =
                        cellVolume * double(inertiaInv[column][row]);
                }
            }
            Eigen::LLT<Matrix6> mobilityFactor(mobility);
            if (mobilityFactor.info() != Eigen::Success) {
                pressureSolveSucceeded = false;
                pressureResidual = std::numeric_limits<float>::infinity();
                return;
            }
            Matrix6 const mobilityRoot = mobilityFactor.matrixL();
            lowRankY = rigidY * mobilityRoot;
            lowRankZ = rigidZ * mobilityRoot;
        }

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> factor;
        factor.compute(baseYY);
        if (factor.info() != Eigen::Success) {
            pressureSolveSucceeded = false;
            pressureResidual = std::numeric_limits<float>::infinity();
            return;
        }

        Eigen::MatrixXd inverseLowRankY;
        Eigen::LDLT<Matrix6> woodburyFactor;
        if (dynamicBody) {
            inverseLowRankY = factor.solve(lowRankY);
            Matrix6 const woodbury =
                Matrix6::Identity()
                + lowRankY.transpose() * inverseLowRankY;
            woodburyFactor.compute(woodbury);
            if (woodburyFactor.info() != Eigen::Success) {
                pressureSolveSucceeded = false;
                pressureResidual = std::numeric_limits<float>::infinity();
                return;
            }
        }

        auto solveUnconstrained = [&](Eigen::MatrixXd const & rhs) {
            Eigen::MatrixXd result = factor.solve(rhs);
            if (dynamicBody) {
                // Woodbury applies the rank-six rigid-body correction exactly.
                result -= inverseLowRankY * woodburyFactor.solve(
                    lowRankY.transpose() * result);
            }
            return result;
        };

        auto multiplyYZ = [&](Eigen::VectorXd const & value) {
            Eigen::VectorXd result =
                Eigen::VectorXd::Zero(unconstrainedCount);
            for (int column = 0; column < candidateCount; ++column) {
                int const row = candidatePressureRows[column];
                if (row >= 0)
                    result[row] +=
                        candidateCoupling[column] * value[column];
            }
            if (dynamicBody)
                result.noalias() +=
                    lowRankY * (lowRankZ.transpose() * value);
            return result;
        };
        auto multiplyYZTranspose = [&](Eigen::VectorXd const & value) {
            Eigen::VectorXd result = Eigen::VectorXd::Zero(candidateCount);
            for (int column = 0; column < candidateCount; ++column) {
                int const row = candidatePressureRows[column];
                if (row >= 0)
                    result[column] =
                        candidateCoupling[column] * value[row];
            }
            if (dynamicBody)
                result.noalias() +=
                    lowRankZ * (lowRankY.transpose() * value);
            return result;
        };
        auto multiplyZZ = [&](Eigen::VectorXd const & value) {
            Eigen::VectorXd result =
                (candidateDiagonal.array() * value.array()).matrix();
            if (dynamicBody)
                result.noalias() +=
                    lowRankZ * (lowRankZ.transpose() * value);
            return result;
        };

        Eigen::VectorXd solvedRhsY = solveUnconstrained(rhsY);

        Eigen::VectorXd candidatePressure =
            Eigen::VectorXd::Zero(candidateCount);
        Eigen::VectorXd reducedGradient =
            Eigen::VectorXd::Zero(candidateCount);
        bool qpConverged = true;
        if (candidateCount > 0) {
            Eigen::VectorXd reducedRhs =
                rhsZ - multiplyYZTranspose(solvedRhsY);

            if (candidateCount <= 64) {
                // Small systems are fastest with the exact dense Schur
                // complement and active-set solve.
                Eigen::MatrixXd denseYZ =
                    Eigen::MatrixXd::Zero(
                        unconstrainedCount,
                        candidateCount);
                for (int column = 0; column < candidateCount; ++column) {
                    int const row = candidatePressureRows[column];
                    if (row >= 0)
                        denseYZ(row, column) = candidateCoupling[column];
                }
                if (dynamicBody)
                    denseYZ.noalias() +=
                        lowRankY * lowRankZ.transpose();
                Eigen::MatrixXd solvedYZ = solveUnconstrained(denseYZ);
                Eigen::MatrixXd denseZZ =
                    candidateDiagonal.asDiagonal();
                if (dynamicBody)
                    denseZZ.noalias() +=
                        lowRankZ * lowRankZ.transpose();
                Eigen::MatrixXd reducedHessian =
                    denseZZ - denseYZ.transpose() * solvedYZ;
                reducedHessian =
                    0.5 * (reducedHessian + reducedHessian.transpose());

                BoundQpResult qpResult =
                    solveNonnegativeQuadratic(reducedHessian, reducedRhs);
                candidatePressure        = std::move(qpResult.Pressure);
                reducedGradient          = std::move(qpResult.Gradient);
                wallSeparationIterations = qpResult.Sweeps;
                qpConverged              = qpResult.Converged;
            } else {
                // Avoid one sparse solve per boundary face. Applying the
                // reduced Hessian matrix-free preserves the same QP:
                // Hq = Azz*q - Azy*Ayy^-1*Ayz*q.
                auto applyReducedHessian =
                    [&](Eigen::VectorXd const & value) {
                        Eigen::VectorXd const eliminated =
                            solveUnconstrained(multiplyYZ(value));
                        Eigen::VectorXd result =
                            multiplyZZ(value)
                            - multiplyYZTranspose(eliminated);
                        return result;
                    };

                Eigen::VectorXd diagonal = candidateDiagonal;
                if (dynamicBody)
                    diagonal.array() +=
                        lowRankZ.rowwise().squaredNorm().array();
                diagonal = diagonal.cwiseMax(1e-8);
                Eigen::VectorXd const scaling = diagonal.cwiseSqrt();
                double lipschitz =
                    (candidateDiagonal.array() / diagonal.array())
                        .maxCoeff();
                if (dynamicBody) {
                    Eigen::MatrixXd scaledLowRank = lowRankZ;
                    scaledLowRank.array().colwise() /= scaling.array();
                    Matrix6 const lowRankGram =
                        scaledLowRank.transpose() * scaledLowRank;
                    Eigen::SelfAdjointEigenSolver<Matrix6> eigenSolver(
                        lowRankGram,
                        Eigen::EigenvaluesOnly);
                    if (eigenSolver.info() == Eigen::Success)
                        lipschitz += eigenSolver.eigenvalues().maxCoeff();
                    else
                        lipschitz += scaledLowRank.squaredNorm();
                }
                lipschitz = std::max(lipschitz, 1e-8);
                Eigen::VectorXd const scaledRhs =
                    (reducedRhs.array() / scaling.array()).matrix();
                auto applyScaledHessian =
                    [&](Eigen::VectorXd const & value) {
                        Eigen::VectorXd const pressure =
                            (value.array() / scaling.array()).matrix();
                        Eigen::VectorXd result =
                            (applyReducedHessian(pressure).array()
                             / scaling.array())
                                .matrix();
                        return result;
                    };

                candidatePressure =
                    (reducedRhs.array() / diagonal.array())
                        .max(0.0)
                        .matrix();
                for (int index = 0; index < candidateCount; ++index) {
                    auto const previous =
                        _wallPressureWarmStart.find(candidateKeys[index]);
                    if (previous != _wallPressureWarmStart.end())
                        candidatePressure[index] =
                            std::max(previous->second, 0.0);
                }
                Eigen::VectorXd scaledPressure =
                    (candidatePressure.array() * scaling.array()).matrix();
                constexpr double MatrixFreeKktTolerance = 1e-6;
                qpConverged = false;
                std::vector<char> freeVariables(candidateCount, 0);
                for (int index = 0; index < candidateCount; ++index) {
                    freeVariables[index] =
                        candidatePressure[index] > MatrixFreeKktTolerance
                        || reducedRhs[index] > MatrixFreeKktTolerance;
                }

                std::unordered_set<std::size_t> visitedSets;
                int hessianApplications = 0;
                for (int activeSetIteration = 0;
                     activeSetIteration < 64;
                     ++activeSetIteration) {
                    std::size_t signature = 1469598103934665603ull;
                    int         freeCount = 0;
                    for (char const isFree : freeVariables) {
                        signature ^= std::size_t(isFree);
                        signature *= 1099511628211ull;
                        freeCount += isFree ? 1 : 0;
                    }
                    if (! visitedSets.insert(signature).second)
                        break;

                    scaledPressure =
                        (candidatePressure.array() * scaling.array()).matrix();
                    Eigen::VectorXd freeRhs = scaledRhs;
                    for (int index = 0; index < candidateCount; ++index) {
                        if (! freeVariables[index]) {
                            scaledPressure[index] = 0.0;
                            freeRhs[index]        = 0.0;
                        }
                    }

                    auto applyFreeHessian =
                        [&](Eigen::VectorXd const & value) {
                            Eigen::VectorXd freeValue = value;
                            for (int index = 0; index < candidateCount; ++index) {
                                if (! freeVariables[index])
                                    freeValue[index] = 0.0;
                            }
                            Eigen::VectorXd result =
                                applyScaledHessian(freeValue);
                            ++hessianApplications;
                            for (int index = 0; index < candidateCount; ++index) {
                                if (! freeVariables[index])
                                    result[index] = 0.0;
                            }
                            return result;
                        };

                    if (freeCount > 0) {
                        Eigen::VectorXd residual =
                            freeRhs - applyFreeHessian(scaledPressure);
                        Eigen::VectorXd direction = residual;
                        double residualSquared = residual.squaredNorm();
                        int const maximumCgIterations =
                            std::min(512, 2 * freeCount + 32);
                        for (int cgIteration = 0;
                             cgIteration < maximumCgIterations
                             && residualSquared > 1e-28;
                             ++cgIteration) {
                            Eigen::VectorXd const hessianDirection =
                                applyFreeHessian(direction);
                            double const curvature =
                                direction.dot(hessianDirection);
                            if (! std::isfinite(curvature)
                                || curvature <= 1e-18)
                                break;

                            double const step = residualSquared / curvature;
                            scaledPressure.noalias() += step * direction;
                            residual.noalias() -= step * hessianDirection;
                            if (! scaledPressure.allFinite()
                                || ! residual.allFinite())
                                break;

                            double const kktResidual =
                                (scaling.array() * residual.array())
                                    .abs()
                                    .maxCoeff();
                            if (kktResidual <= MatrixFreeKktTolerance)
                                break;

                            double const nextResidualSquared =
                                residual.squaredNorm();
                            direction =
                                residual
                                + (nextResidualSquared / residualSquared)
                                    * direction;
                            residualSquared = nextResidualSquared;
                        }
                    } else {
                        scaledPressure.setZero();
                    }

                    Eigen::VectorXd trialPressure =
                        (scaledPressure.array() / scaling.array()).matrix();
                    bool removedNegativePressure = false;
                    for (int index = 0; index < candidateCount; ++index) {
                        if (freeVariables[index]
                            && trialPressure[index]
                                < -MatrixFreeKktTolerance) {
                            freeVariables[index] = 0;
                            candidatePressure[index] = 0.0;
                            removedNegativePressure = true;
                        } else if (freeVariables[index]) {
                            candidatePressure[index] =
                                std::max(trialPressure[index], 0.0);
                        } else {
                            candidatePressure[index] = 0.0;
                        }
                    }
                    if (removedNegativePressure)
                        continue;

                    reducedGradient =
                        applyReducedHessian(candidatePressure) - reducedRhs;
                    ++hessianApplications;
                    if (! reducedGradient.allFinite())
                        break;

                    bool   releasedConstraint = false;
                    double kktResidual        = 0.0;
                    for (int index = 0; index < candidateCount; ++index) {
                        if (! freeVariables[index]
                            && reducedGradient[index]
                                < -MatrixFreeKktTolerance) {
                            freeVariables[index] = 1;
                            releasedConstraint   = true;
                        }
                        double const violation =
                            candidatePressure[index]
                                    > MatrixFreeKktTolerance
                            ? std::abs(reducedGradient[index])
                            : std::max(-reducedGradient[index], 0.0);
                        kktResidual = std::max(kktResidual, violation);
                    }
                    if (! releasedConstraint
                        && kktResidual <= MatrixFreeKktTolerance) {
                        qpConverged = true;
                        break;
                    }
                }
                if (! qpConverged) {
                    scaledPressure =
                        (candidatePressure.array() * scaling.array()).matrix();
                    Eigen::VectorXd accelerated = scaledPressure;
                    double          momentum     = 1.0;
                    int const maximumIterations =
                        std::min(2000, 8 * candidateCount + 128);
                    for (int iteration = 0;
                         iteration < maximumIterations;
                         ++iteration) {
                        Eigen::VectorXd const acceleratedGradient =
                            applyScaledHessian(accelerated) - scaledRhs;
                        ++hessianApplications;
                        Eigen::VectorXd const nextScaledPressure =
                            (accelerated.array()
                             - acceleratedGradient.array() / lipschitz)
                                .max(0.0)
                                .matrix();
                        if (! acceleratedGradient.allFinite()
                            || ! nextScaledPressure.allFinite())
                            break;

                        if ((iteration + 1) % 4 == 0
                            || iteration + 1 == maximumIterations) {
                            candidatePressure =
                                (nextScaledPressure.array() / scaling.array())
                                    .matrix();
                            reducedGradient =
                                applyReducedHessian(candidatePressure)
                                - reducedRhs;
                            ++hessianApplications;
                            if (! reducedGradient.allFinite())
                                break;
                            double kktResidual = 0.0;
                            for (int index = 0; index < candidateCount; ++index) {
                                double const violation =
                                    candidatePressure[index]
                                            > MatrixFreeKktTolerance
                                    ? std::abs(reducedGradient[index])
                                    : std::max(-reducedGradient[index], 0.0);
                                kktResidual =
                                    std::max(kktResidual, violation);
                            }
                            if (kktResidual <= MatrixFreeKktTolerance) {
                                qpConverged = true;
                                break;
                            }
                        }

                        double const nextMomentum =
                            0.5 * (1.0 + std::sqrt(
                                1.0 + 4.0 * momentum * momentum));
                        Eigen::VectorXd nextAccelerated =
                            nextScaledPressure
                            + ((momentum - 1.0) / nextMomentum)
                                * (nextScaledPressure - scaledPressure);
                        if ((accelerated - nextScaledPressure).dot(
                                nextScaledPressure - scaledPressure)
                            > 0.0) {
                            nextAccelerated = nextScaledPressure;
                            momentum        = 1.0;
                        } else {
                            momentum = nextMomentum;
                        }
                        scaledPressure = nextScaledPressure;
                        accelerated    = std::move(nextAccelerated);
                    }

                    if (! qpConverged) {
                        candidatePressure =
                            (scaledPressure.array() / scaling.array()).matrix();
                        reducedGradient =
                            applyReducedHessian(candidatePressure) - reducedRhs;
                        ++hessianApplications;
                    }
                }
                wallSeparationIterations = hessianApplications;
            }
        } else {
            wallSeparationIterations = 1;
        }

        _wallPressureWarmStart.clear();
        _wallPressureWarmStart.reserve(candidateCount);
        for (int index = 0; index < candidateCount; ++index)
            _wallPressureWarmStart.emplace(
                candidateKeys[index],
                candidatePressure[index]);

        Eigen::VectorXd unconstrainedPressure = solveUnconstrained(
            rhsY - multiplyYZ(candidatePressure));
        if (! unconstrainedPressure.allFinite()
            || ! candidatePressure.allFinite()) {
            pressureSolveSucceeded = false;
            pressureResidual = std::numeric_limits<float>::infinity();
            return;
        }

        Eigen::VectorXd unconstrainedResidual =
            baseYY * unconstrainedPressure
            + multiplyYZ(candidatePressure)
            - rhsY;
        if (dynamicBody) {
            unconstrainedResidual +=
                lowRankY
                * (lowRankY.transpose() * unconstrainedPressure);
        }
        pressureResidual = float(
            unconstrainedResidual.norm()
            / std::max(rhsY.norm(), 1.0));

        Eigen::VectorXd pressure =
            unconstrainedPressure.head(pressureCount);
        for (int row = 0; row < pressureCount; ++row)
            m_p[rowToCell[row]] = float(pressure[row]);

        auto boundaryPressure = [&](WallFace const & wallFace) {
            if (wallFace.Candidate)
                return float(std::max(
                    candidatePressure[wallFace.CandidateIndex],
                    0.0));
            return float(
                unconstrainedPressure[wallFace.UnconstrainedIndex]);
        };

        m_vel = intermediateVelocity;
        for (int row = 0; row < pressureCount; ++row) {
            glm::ivec3 const cell = decodeCell(rowToCell[row]);
            float const pressureValue = float(pressure[row]);
            for (FaceNeighbor const & side : FaceNeighbors) {
                glm::ivec3 const face     = cell + side.FaceOffset;
                glm::ivec3 const neighbor = cell + side.CellOffset;
                int const        faceIdx  = gridOffset(face);
                float const openWeight =
                    faceWeight(side.Direction, faceIdx);
                if (openWeight <= 1e-6f
                    || ! isValidCell(neighbor)
                    || m_s[gridOffset(neighbor)] <= 0.0f)
                    continue;

                float pressureScale = 1.0f;
                int const neighborIdx = gridOffset(neighbor);
                if (cellToRow[neighborIdx] < 0
                    && m_type[neighborIdx] == EMPTY_CELL) {
                    pressureScale = ghostFluidPressureScale(
                        liquidLevelSet,
                        rowToCell[row],
                        neighborIdx);
                }
                m_vel[faceIdx][side.Direction] +=
                    side.DivergenceSign
                    * pressureScale
                    * pressureValue;
            }
        }

        std::vector<BoundaryPressureSample> boundarySamples;
        boundarySamples.reserve(wallFaces.size());
        for (WallFace const & wallFace : wallFaces) {
            float const q = boundaryPressure(wallFace);
            if (wallFace.OpenWeight <= 1e-6f) {
                m_vel[wallFace.FaceIndex][wallFace.Direction] =
                    intermediateVelocity[wallFace.FaceIndex][wallFace.Direction]
                    + wallFace.DivergenceSign
                        * wallFace.GhostScale
                        * (float(pressure[wallFace.Row]) - q);
            }
            if (wallFace.BodyBoundary) {
                boundarySamples.push_back(BoundaryPressureSample {
                    .FaceIndex        = wallFace.FaceIndex,
                    .Direction        = wallFace.Direction,
                    .DivergenceSign   = wallFace.DivergenceSign,
                    .BoundaryWeight   = wallFace.BoundaryWeight,
                    .BoundaryPressure = q,
                    .Contact          = ! wallFace.Candidate || q > 1e-7f,
                });
            }
        }
        applyRigidBodyFeedback(boundarySamples, dt);

        glm::vec3 predictedLinearDelta(0.0f);
        glm::vec3 predictedAngularDelta(0.0f);
        if (dynamicBody) {
            predictedLinearDelta =
                m_feedbackForce / m_body->mass * dt;
            predictedAngularDelta =
                m_body->GetInertiaWorldInv() * m_feedbackTorque * dt;
        }

        wallSeparationKktResidual   = 0.0f;
        wallSeparationQpResidual    = 0.0f;
        wallSeparationSpeedResidual = 0.0f;
        maximumBoundaryPressure     = 0.0f;
        for (WallFace const & wallFace : wallFaces) {
            if (! wallFace.Candidate)
                continue;

            float const q = boundaryPressure(wallFace);
            maximumBoundaryPressure =
                std::max(maximumBoundaryPressure, q);
            float const boundaryVelocity =
                intermediateVelocity[wallFace.FaceIndex][wallFace.Direction]
                + wallFace.DivergenceSign
                    * wallFace.GhostScale
                    * (float(pressure[wallFace.Row]) - q);
            float finalWallVelocity = wallFace.WallVelocity;
            if (dynamicBody && wallFace.BodyBoundary) {
                glm::ivec3 const face = decodeCell(wallFace.FaceIndex);
                glm::vec3 const lever =
                    faceCenter(face, wallFace.Direction) - m_body->position;
                finalWallVelocity += (
                    predictedLinearDelta
                    + glm::cross(predictedAngularDelta, lever))
                                         [wallFace.Direction];
            }
            float const separationSpeed =
                -wallFace.DivergenceSign
                * (boundaryVelocity - finalWallVelocity);
            float violation = std::max(-q, 0.0f);
            if (q > 1e-6f)
                violation = std::max(violation, std::abs(separationSpeed));
            else
                violation = std::max(
                    violation,
                    std::max(-separationSpeed, 0.0f));
            violation = std::max(
                violation,
                std::abs(q * separationSpeed));
            wallSeparationSpeedResidual =
                std::max(wallSeparationSpeedResidual, violation);
        }
        for (int index = 0; index < candidateCount; ++index) {
            float violation =
                float(std::max(-candidatePressure[index], 0.0));
            if (candidatePressure[index] > 1e-7)
                violation = std::max(
                    violation,
                    float(std::abs(reducedGradient[index])));
            else
                violation = std::max(
                    violation,
                    float(std::max(-reducedGradient[index], 0.0)));
            wallSeparationQpResidual =
                std::max(wallSeparationQpResidual, violation);
        }
        wallSeparationKktResidual = std::max(
            wallSeparationSpeedResidual,
            wallSeparationQpResidual);

        pressureSolveSucceeded =
            std::isfinite(pressureResidual)
            && pressureResidual < 2e-4f
            && wallSeparationKktResidual < 2e-4f;
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

            glm::ivec3 const face {
                sample.FaceIndex % m_iCellX,
                (sample.FaceIndex / m_iCellX) % m_iCellY,
                sample.FaceIndex / (m_iCellX * m_iCellY),
            };
            glm::vec3 const applicationPoint =
                faceCenter(face, sample.Direction);
            glm::vec3 forceDirection(0.0f);
            forceDirection[sample.Direction] = sample.DivergenceSign;
            float const impulseScale =
                sample.BoundaryWeight * m_h * m_h * m_h;
            glm::vec3 const force =
                forceDirection
                * sample.BoundaryPressure
                * impulseScale
                * invDt;
            m_feedbackForce += force;
            m_feedbackTorque +=
                glm::cross(applicationPoint - m_body->position, force);
        }
    }
} // namespace VCX::Labs::Final
