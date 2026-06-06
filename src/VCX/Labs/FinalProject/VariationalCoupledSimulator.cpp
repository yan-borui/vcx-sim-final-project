#include "Labs/FinalProject/VariationalCoupledSimulator.h"

#include <algorithm>
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

    void VariationalCoupledSimulator::setupScene(int res) {
        Simulator::setupScene(res);
        for (auto & fractions : _faceFluidFraction)
            fractions.assign(m_iNumCells, 0.0f);
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
        bool touchesAir = false;

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
        int dir) const {
        if (! isTankFaceOpen(face, dir))
            return 0.0f;
        if (! m_body)
            return 1.0f;

        glm::vec3 const center = faceCenter(face, dir);
        float const centerDistance = m_body->GetSDF(center);
        float const sampleRadius = 0.5f * std::sqrt(3.0f) * m_h;
        if (centerDistance >= sampleRadius)
            return 1.0f;
        if (centerDistance <= -sampleRadius)
            return 0.0f;

        int const sampleCount = std::max(volumeSamplesPerAxis, 1);
        int fluidSamples = 0;

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
                    int const idx = gridOffset(face);
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
        int dir,
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
        int numIters,
        float dt,
        float overRelaxation,
        bool compensateDrift) {
        (void) dt;
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
        m_feedbackForce = glm::vec3(0.0f);
        m_feedbackTorque = glm::vec3(0.0f);
        if (rowToCell.empty())
            return;

        auto decodeCell = [&](int idx) {
            return glm::ivec3(
                idx % m_iCellX,
                (idx / m_iCellX) % m_iCellY,
                idx / (m_iCellX * m_iCellY));
        };

        std::vector<char> candidate(rowToCell.size(), 0);
        std::vector<char> separating(rowToCell.size(), 0);
        for (int row = 0; row < int(rowToCell.size()); ++row)
            candidate[row] = isWallSeparationCandidate(decodeCell(rowToCell[row])) ? 1 : 0;

        Eigen::VectorXd pressure = Eigen::VectorXd::Zero(int(rowToCell.size()));
        auto solveWithActiveSet = [&]() {
            int const matrixSize = int(rowToCell.size());
            Eigen::SparseMatrix<double> matrix(matrixSize, matrixSize);
            std::vector<Eigen::Triplet<double>> triplets;
            triplets.reserve(matrixSize * 7);
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(matrixSize);

            for (int row = 0; row < matrixSize; ++row) {
                if (separating[row]) {
                    triplets.emplace_back(row, row, 1.0);
                    continue;
                }

                int const idx = rowToCell[row];
                glm::ivec3 const cell = decodeCell(idx);
                double diagonal = 0.0;
                double weightedDivergence = 0.0;

                for (FaceNeighbor const & side : FaceNeighbors) {
                    glm::ivec3 const face = cell + side.FaceOffset;
                    glm::ivec3 const neighbor = cell + side.CellOffset;
                    int const faceIdx = gridOffset(face);
                    float const openWeight = faceWeight(side.Direction, faceIdx);

                    if (isSolidPressureCell(neighbor)) {
                        float const contactWeight = 1.0f - openWeight;
                        if (contactWeight <= 1e-6f)
                            continue;
                        float const relativeVelocity =
                            m_vel[faceIdx][side.Direction]
                            - solidVelocity(face, side.Direction, neighbor);
                        weightedDivergence += side.DivergenceSign
                            * double(contactWeight) * double(relativeVelocity);
                        diagonal += contactWeight;
                        continue;
                    }

                    if (openWeight <= 1e-6f)
                        continue;

                    weightedDivergence += side.DivergenceSign
                        * double(openWeight) * double(m_vel[faceIdx][side.Direction]);
                    diagonal += openWeight;

                    int const neighborRow = cellToRow[gridOffset(neighbor)];
                    if (neighborRow >= 0 && ! separating[neighborRow])
                        triplets.emplace_back(row, neighborRow, -double(openWeight));
                }

                if (compensateDrift && m_particleRestDensity > 0.0f) {
                    float const densityError = m_particleDensity[idx] - m_particleRestDensity;
                    if (densityError > 0.0f)
                        weightedDivergence -= compensateDriftWeight * densityError;
                }

                if (diagonal <= 1e-8) {
                    diagonal = 1.0;
                    weightedDivergence = 0.0;
                }
                triplets.emplace_back(row, row, diagonal);
                rhs[row] = -weightedDivergence;
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

            for (int row = 0; row < int(separating.size()); ++row) {
                if (separating[row])
                    pressure[row] = 0.0;
            }
        };

        bool changedActiveSet = false;
        int const maxActiveSetIterations = enableWallSeparation
            ? std::min<int>(32, int(rowToCell.size()) + 1)
            : 1;
        for (int iter = 0; iter < maxActiveSetIterations; ++iter) {
            solveWithActiveSet();
            changedActiveSet = false;
            if (! enableWallSeparation)
                break;

            for (int row = 0; row < int(rowToCell.size()); ++row) {
                if (candidate[row] && ! separating[row] && pressure[row] < -1e-5) {
                    separating[row] = 1;
                    changedActiveSet = true;
                }
            }
            if (! changedActiveSet)
                break;
        }
        if (changedActiveSet)
            solveWithActiveSet();

        for (int row = 0; row < int(rowToCell.size()); ++row)
            m_p[rowToCell[row]] = float(pressure[row]);

        for (int row = 0; row < int(rowToCell.size()); ++row) {
            glm::ivec3 const cell = decodeCell(rowToCell[row]);
            float const pressureValue = float(pressure[row]);

            for (FaceNeighbor const & side : FaceNeighbors) {
                glm::ivec3 const face = cell + side.FaceOffset;
                glm::ivec3 const neighbor = cell + side.CellOffset;
                int const faceIdx = gridOffset(face);
                float const weight = isSolidPressureCell(neighbor)
                    ? 1.0f - faceWeight(side.Direction, faceIdx)
                    : faceWeight(side.Direction, faceIdx);
                if (weight > 1e-6f)
                    m_vel[faceIdx][side.Direction] += side.DivergenceSign * pressureValue;
            }
        }

        enforceSolidBoundaryVelocities(rowToCell, separating);
        applyRigidBodyFeedback(rowToCell);
    }

    void VariationalCoupledSimulator::enforceSolidBoundaryVelocities(
        std::vector<int> const & pressureCells,
        std::vector<char> const & separatingCells) {
        for (int row = 0; row < int(pressureCells.size()); ++row) {
            int const idx = pressureCells[row];
            glm::ivec3 const cell {
                idx % m_iCellX,
                (idx / m_iCellX) % m_iCellY,
                idx / (m_iCellX * m_iCellY),
            };

            for (FaceNeighbor const & side : FaceNeighbors) {
                glm::ivec3 const neighbor = cell + side.CellOffset;
                if (! isSolidPressureCell(neighbor))
                    continue;

                glm::ivec3 const face = cell + side.FaceOffset;
                int const faceIdx = gridOffset(face);
                if (1.0f - faceWeight(side.Direction, faceIdx) <= 1e-6f)
                    continue;

                float const wallVelocity = solidVelocity(face, side.Direction, neighbor);
                float & velocity = m_vel[faceIdx][side.Direction];
                if (! enableWallSeparation || ! separatingCells[row]) {
                    velocity = wallVelocity;
                } else if (side.DivergenceSign < 0.0f) {
                    velocity = std::max(velocity, wallVelocity);
                } else {
                    velocity = std::min(velocity, wallVelocity);
                }
            }
        }
    }

    void VariationalCoupledSimulator::applyRigidBodyFeedback(
        std::vector<int> const & pressureCells) {
        if (! m_body)
            return;

        for (int idx : pressureCells) {
            float const pressure = m_p[idx];
            if (pressure <= 0.0f)
                continue;

            glm::ivec3 const cell {
                idx % m_iCellX,
                (idx / m_iCellX) % m_iCellY,
                idx / (m_iCellX * m_iCellY),
            };
            glm::vec3 const cellPosition = pressureCellCenter(cell);
            if (m_body->GetSDF(cellPosition) >= m_h)
                continue;

            float const epsilon = 0.001f;
            glm::vec3 const gradient {
                m_body->GetSDF(cellPosition + glm::vec3(epsilon, 0, 0))
                    - m_body->GetSDF(cellPosition - glm::vec3(epsilon, 0, 0)),
                m_body->GetSDF(cellPosition + glm::vec3(0, epsilon, 0))
                    - m_body->GetSDF(cellPosition - glm::vec3(0, epsilon, 0)),
                m_body->GetSDF(cellPosition + glm::vec3(0, 0, epsilon))
                    - m_body->GetSDF(cellPosition - glm::vec3(0, 0, epsilon)),
            };
            float const gradientLength = glm::length(gradient);
            if (gradientLength <= 1e-6f)
                continue;

            glm::vec3 const normal = gradient / gradientLength;
            float const area = m_h * m_h;
            glm::vec3 const force = -normal * pressure * area * 45.0f;
            m_feedbackForce += force;
            m_feedbackTorque += glm::cross(cellPosition - m_body->position, force);
        }
    }
}
