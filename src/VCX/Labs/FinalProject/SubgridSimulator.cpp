#include "Labs/FinalProject/SubgridSimulator.h"

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

    void SubgridSimulator::setupScene(int res) {
        Simulator::setupScene(res);

        voxelizeDynamicBody = false;
        gravity           = { 6.0f, 0.0f, 0.0f };
        m_fRatio          = 0.9f;
        numPressureIters  = 120;
        separateParticles = false;
        compensateDrift   = false;

        setupChannelParticles();

        for (auto & fractions : _faceOpenFraction)
            fractions.assign(m_iNumCells, 0.0f);
        for (auto & fractions : _faceFluidFraction)
            fractions.assign(m_iNumCells, 0.0f);
        pressureResidual       = 0.0f;
        pressureSolveSucceeded = true;
        updateFaceFluidFractions(buildParticleLevelSet());
    }

    void SubgridSimulator::setupChannelParticles() {
        float const wallMin    = -0.5f + m_h;
        float const wallMax    = 0.5f - m_h;
        float const spacingX   = 0.50f * m_h;
        float const spacingY   = 0.46f * m_h;
        float const spacingZ   = 0.50f * m_h;
        float const sourceMaxX = m_body
            ? m_body->position.x - 0.5f * m_body->dim.x - 0.75f * m_h
            : -0.05f;

        std::vector<glm::vec3> positions;
        for (float x = wallMin + m_particleRadius; x <= sourceMaxX; x += spacingX) {
            int const xLayer = int(std::round((x - wallMin) / spacingX));
            for (float y = wallMin + m_particleRadius; y <= wallMax - m_particleRadius; y += spacingY) {
                int const   yLayer  = int(std::round((y - wallMin) / spacingY));
                float const zOffset = ((xLayer + yLayer) & 1) ? 0.25f * spacingZ : 0.0f;
                for (float z = wallMin + m_particleRadius + zOffset; z <= wallMax - m_particleRadius; z += spacingZ) {
                    glm::vec3 const position { x, y, z };
                    if (! m_body || m_body->GetSDF(position) > m_particleRadius)
                        positions.push_back(position);
                }
            }
        }

        m_iNumSpheres = int(positions.size());
        m_particlePos = std::move(positions);
        m_particleVel.assign(m_iNumSpheres, glm::vec3(0.0f));
        m_particleColor.assign(m_iNumSpheres, glm::vec3(1.0f));
        m_hashtable.assign(m_iNumSpheres, 0);
        m_hashtableindex.assign(m_iNumCells + 1, 0);
        m_particleRestDensity = 0.0f;
    }

    bool SubgridSimulator::isValidCell(glm::ivec3 const & cell) const {
        return cell.x >= 0 && cell.x < m_iCellX
            && cell.y >= 0 && cell.y < m_iCellY
            && cell.z >= 0 && cell.z < m_iCellZ;
    }

    glm::vec3 SubgridSimulator::pressureCellCenter(glm::ivec3 const & cell) const {
        return (glm::vec3(cell) + glm::vec3(0.5f)) * m_h - glm::vec3(0.5f);
    }

    bool SubgridSimulator::isTankFaceOpen(glm::ivec3 const & face, int dir) {
        glm::ivec3 left = face;
        left[dir] -= 1;
        glm::ivec3 const right = face;
        if (! isValidCell(left) || ! isValidCell(right))
            return false;
        return m_s[index2GridOffset(left)] > 0.0f
            && m_s[index2GridOffset(right)] > 0.0f;
    }

    bool SubgridSimulator::isTankSolidCell(glm::ivec3 const & cell) {
        if (! isValidCell(cell) || m_s[index2GridOffset(cell)] <= 0.0f)
            return true;
        return false;
    }

    bool SubgridSimulator::isBodyInteriorPressureCell(glm::ivec3 const & cell) {
        return ! isTankSolidCell(cell)
            && m_body
            && m_body->GetSDF(pressureCellCenter(cell)) < -1e-5f;
    }

    bool SubgridSimulator::hasFluidSupportAcrossOpenFace(glm::ivec3 const & cell) {
        if (! isBodyInteriorPressureCell(cell))
            return false;

        for (FaceNeighbor const & side : FaceNeighbors) {
            glm::ivec3 const face = cell + side.FaceOffset;
            if (! isValidCell(face))
                continue;

            int const faceIdx = index2GridOffset(face);
            if (_faceFluidFraction[side.Direction][faceIdx] <= 1e-6f)
                continue;

            glm::ivec3 const neighbor = cell + side.CellOffset;
            if (isTankSolidCell(neighbor))
                continue;

            int const neighborIdx = index2GridOffset(neighbor);
            if (m_type[neighborIdx] == FLUID_CELL)
                return true;
        }

        return false;
    }

    bool SubgridSimulator::isPressureUnknownCell(glm::ivec3 const & cell) {
        if (isTankSolidCell(cell))
            return false;

        int const idx = index2GridOffset(cell);
        return m_type[idx] == FLUID_CELL
            || hasFluidSupportAcrossOpenFace(cell);
    }

    glm::vec2 SubgridSimulator::estimateFaceFractions(
        glm::ivec3 const &         face,
        int                        dir,
        std::vector<float> const & liquidLevelSet) {
        if (! isTankFaceOpen(face, dir))
            return glm::vec2(0.0f);

        glm::vec3 offset(0.5f);
        offset[dir]                  = 0.0f;
        glm::vec3 const center       = (glm::vec3(face) + offset) * m_h - glm::vec3(0.5f);
        int const       sampleCount  = std::clamp(volumeSamplesPerAxis, 1, 16);
        int             openSamples  = 0;
        int             fluidSamples = 0;

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

    void SubgridSimulator::updateFaceFluidFractions(
        std::vector<float> const & liquidLevelSet) {
        partiallyOpenFaceCount   = 0;
        partiallyFilledFaceCount = 0;
        minimumOpenFaceFraction  = 1.0f;
        minimumFluidMassFraction = 1.0f;

        for (int k = 0; k < m_iCellZ; ++k) {
            for (int j = 0; j < m_iCellY; ++j) {
                for (int i = 0; i < m_iCellX; ++i) {
                    glm::ivec3 const face { i, j, k };
                    int const        idx = index2GridOffset(face);
                    for (int dir = 0; dir < 3; ++dir) {
                        glm::vec2 const fractions =
                            estimateFaceFractions(face, dir, liquidLevelSet);
                        _faceOpenFraction[dir][idx]  = fractions.x;
                        _faceFluidFraction[dir][idx] = fractions.y;
                        if (fractions.x > 0.0f && fractions.x < 1.0f) {
                            ++partiallyOpenFaceCount;
                            minimumOpenFaceFraction =
                                std::min(minimumOpenFaceFraction, fractions.x);
                        }
                        if (fractions.y > 0.0f && fractions.y < 1.0f) {
                            ++partiallyFilledFaceCount;
                            minimumFluidMassFraction =
                                std::min(minimumFluidMassFraction, fractions.y);
                        }
                    }
                }
            }
        }

        if (partiallyOpenFaceCount == 0)
            minimumOpenFaceFraction = 0.0f;
        if (partiallyFilledFaceCount == 0)
            minimumFluidMassFraction = 0.0f;
    }

    float SubgridSimulator::faceWeight(int dir, int idx) const {
        float const fraction = _faceFluidFraction[dir][idx];
        if (useSubgridWeights)
            return fraction;
        return fraction > 0.5f ? 1.0f : 0.0f;
    }

    void SubgridSimulator::solveIncompressibility(
        int   numIters,
        float dt,
        float overRelaxation,
        bool  compensateDrift) {
        (void) dt;
        (void) overRelaxation;
        (void) compensateDrift;

        std::vector<float> const liquidLevelSet = buildParticleLevelSet();
        updateFaceFluidFractions(liquidLevelSet);

        std::vector<int> cellToMatrix(m_iNumCells, -1);
        std::vector<int> matrixToCell;
        matrixToCell.reserve(m_iNumCells / 3);

        for (int idx = 0; idx < m_iNumCells; ++idx) {
            glm::ivec3 const cell {
                idx % m_iCellX,
                (idx / m_iCellX) % m_iCellY,
                idx / (m_iCellX * m_iCellY),
            };
            if (isPressureUnknownCell(cell)) {
                cellToMatrix[idx] = int(matrixToCell.size());
                matrixToCell.push_back(idx);
            }
        }

        int const matrixSize = int(matrixToCell.size());
        std::fill(m_p.begin(), m_p.end(), 0.0f);
        m_feedbackForce        = glm::vec3(0.0f);
        m_feedbackTorque       = glm::vec3(0.0f);
        pressureResidual       = 0.0f;
        pressureSolveSucceeded = true;
        if (matrixSize == 0)
            return;

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

                int const        idx = matrixToCell[row];
                glm::ivec3 const cell {
                    idx % m_iCellX,
                    (idx / m_iCellX) % m_iCellY,
                    idx / (m_iCellX * m_iCellY),
                };
                for (FaceNeighbor const & side : FaceNeighbors) {
                    glm::ivec3 const face     = cell + side.FaceOffset;
                    glm::ivec3 const neighbor = cell + side.CellOffset;
                    int const        faceIdx  = index2GridOffset(face);
                    if (faceWeight(side.Direction, faceIdx) <= 1e-6f
                        || isTankSolidCell(neighbor))
                        continue;

                    int const neighborColumn =
                        cellToMatrix[index2GridOffset(neighbor)];
                    if (neighborColumn >= 0) {
                        if (! visitedRows[neighborColumn]) {
                            visitedRows[neighborColumn] = 1;
                            pending.push_back(neighborColumn);
                        }
                    } else if (m_type[index2GridOffset(neighbor)] == EMPTY_CELL) {
                        hasDirichletBoundary = true;
                    }
                }
            }

            if (! hasDirichletBoundary)
                pinnedRows[component.front()] = 1;
        }

        Eigen::SparseMatrix<double>         matrix(matrixSize, matrixSize);
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(matrixSize * 7);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(matrixSize);

        for (int row = 0; row < matrixSize; ++row) {
            if (pinnedRows[row]) {
                triplets.emplace_back(row, row, 1.0);
                continue;
            }

            int const        idx = matrixToCell[row];
            glm::ivec3 const cell {
                idx % m_iCellX,
                (idx / m_iCellX) % m_iCellY,
                idx / (m_iCellX * m_iCellY),
            };

            double diagonal           = 0.0;
            double weightedDivergence = 0.0;
            for (FaceNeighbor const & side : FaceNeighbors) {
                glm::ivec3 const face     = cell + side.FaceOffset;
                int const        faceIdx  = index2GridOffset(face);
                float const      weight   = faceWeight(side.Direction, faceIdx);
                glm::ivec3 const neighbor = cell + side.CellOffset;

                if (weight <= 1e-6f || isTankSolidCell(neighbor)) {
                    m_vel[faceIdx][side.Direction] = 0.0f;
                    continue;
                }

                weightedDivergence += side.DivergenceSign
                    * double(weight) * double(m_vel[faceIdx][side.Direction]);

                int const neighborIdx    = index2GridOffset(neighbor);
                int const neighborColumn = cellToMatrix[neighborIdx];
                if (neighborColumn >= 0) {
                    diagonal += weight;
                    if (! pinnedRows[neighborColumn])
                        triplets.emplace_back(row, neighborColumn, -double(weight));
                } else if (m_type[neighborIdx] == EMPTY_CELL) {
                    diagonal += weight * ghostFluidPressureScale(
                        liquidLevelSet,
                        idx,
                        neighborIdx);
                } else {
                    diagonal += weight;
                }
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
        Eigen::VectorXd pressure = Eigen::VectorXd::Zero(matrixSize);
        if (solver.info() == Eigen::Success)
            pressure = solver.solve(rhs);
        pressureSolveSucceeded =
            solver.info() == Eigen::Success && pressure.allFinite();
        if (! pressureSolveSucceeded) {
            pressureResidual = std::numeric_limits<float>::infinity();
            return;
        }
        pressureResidual = float(solver.error());

        for (int row = 0; row < matrixSize; ++row)
            m_p[matrixToCell[row]] = float(pressure[row]);

        // Each pressure contributes to its six surrounding MAC faces. Adjacent
        // pressure contributions combine into the usual finite-difference gradient.
        for (int row = 0; row < matrixSize; ++row) {
            int const        idx = matrixToCell[row];
            glm::ivec3 const cell {
                idx % m_iCellX,
                (idx / m_iCellX) % m_iCellY,
                idx / (m_iCellX * m_iCellY),
            };
            float const pressureValue = m_p[idx];

            for (FaceNeighbor const & side : FaceNeighbors) {
                glm::ivec3 const face     = cell + side.FaceOffset;
                int const        faceIdx  = index2GridOffset(face);
                glm::ivec3 const neighbor = cell + side.CellOffset;
                if (faceWeight(side.Direction, faceIdx) <= 1e-6f || isTankSolidCell(neighbor))
                    continue;

                int const   neighborIdx = index2GridOffset(neighbor);
                float const pressureScale =
                    cellToMatrix[neighborIdx] < 0 && m_type[neighborIdx] == EMPTY_CELL
                    ? ghostFluidPressureScale(liquidLevelSet, idx, neighborIdx)
                    : 1.0f;
                m_vel[faceIdx][side.Direction] +=
                    side.DivergenceSign * pressureScale * pressureValue;
            }
        }
    }
} // namespace VCX::Labs::Final
