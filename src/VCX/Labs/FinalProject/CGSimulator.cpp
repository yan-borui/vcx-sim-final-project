#include "Labs/FinalProject/CGSimulator.h"

#include <cmath>
#include <limits>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

namespace VCX::Labs::Final {
    namespace {
        struct Neighbor {
            glm::ivec3 CellOffset;
            glm::ivec3 FaceOffset;
            int        Direction;
            float      Sign;
        };

        constexpr Neighbor Neighbors[] = {
            { { -1, 0, 0 }, { 0, 0, 0 }, 0, -1.0f },
            { {  1, 0, 0 }, { 1, 0, 0 }, 0,  1.0f },
            { { 0, -1, 0 }, { 0, 0, 0 }, 1, -1.0f },
            { {  0, 1, 0 }, { 0, 1, 0 }, 1,  1.0f },
            { { 0, 0, -1 }, { 0, 0, 0 }, 2, -1.0f },
            { {  0, 0, 1 }, { 0, 0, 1 }, 2,  1.0f },
        };

        bool finite(glm::vec3 const & v) {
            return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
        }
    }

    void CGSimulator::solveIncompressibility(
        int   numIters,
        float dt,
        float overRelaxation,
        bool  compensateDrift) {
        (void) dt;
        (void) overRelaxation;

        std::vector<int> cellToRow(m_iNumCells, -1);
        std::vector<int> rowToCell;
        rowToCell.reserve(m_iNumCells / 3);

        for (int idx = 0; idx < m_iNumCells; ++idx) {
            if (m_type[idx] != FLUID_CELL)
                continue;
            cellToRow[idx] = int(rowToCell.size());
            rowToCell.push_back(idx);
        }

        std::fill(m_p.begin(), m_p.end(), 0.0f);
        m_feedbackForce  = glm::vec3(0.0f);
        m_feedbackTorque = glm::vec3(0.0f);
        if (rowToCell.empty())
            return;

        auto decode = [&](int idx) {
            return glm::ivec3(
                idx % m_iCellX,
                (idx / m_iCellX) % m_iCellY,
                idx / (m_iCellX * m_iCellY));
        };

        int const matrixSize = int(rowToCell.size());
        std::vector<char> pinned(matrixSize, 0);
        std::vector<char> visited(matrixSize, 0);
        for (int seed = 0; seed < matrixSize; ++seed) {
            if (visited[seed])
                continue;

            bool             hasAirBoundary = false;
            std::vector<int> component;
            std::vector<int> stack { seed };
            visited[seed] = 1;

            while (! stack.empty()) {
                int const row = stack.back();
                stack.pop_back();
                component.push_back(row);

                glm::ivec3 const cell = decode(rowToCell[row]);
                for (Neighbor const & side : Neighbors) {
                    glm::ivec3 const neighbor = cell + side.CellOffset;
                    if (! isValidCell(neighbor))
                        continue;

                    int const neighborIdx = index2GridOffset(neighbor);
                    if (m_s[neighborIdx] <= 0.0f)
                        continue;

                    int const neighborRow = cellToRow[neighborIdx];
                    if (neighborRow >= 0) {
                        if (! visited[neighborRow]) {
                            visited[neighborRow] = 1;
                            stack.push_back(neighborRow);
                        }
                    } else if (m_type[neighborIdx] == EMPTY_CELL) {
                        hasAirBoundary = true;
                    }
                }
            }

            if (! hasAirBoundary && ! component.empty())
                pinned[component.front()] = 1;
        }

        Eigen::SparseMatrix<double>         matrix(matrixSize, matrixSize);
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(matrixSize * 7);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(matrixSize);

        for (int row = 0; row < matrixSize; ++row) {
            if (pinned[row]) {
                triplets.emplace_back(row, row, 1.0);
                continue;
            }

            int const        idx  = rowToCell[row];
            glm::ivec3 const cell = decode(idx);

            double diagonal   = 0.0;
            double divergence = 0.0;

            for (Neighbor const & side : Neighbors) {
                glm::ivec3 const face     = cell + side.FaceOffset;
                glm::ivec3 const neighbor = cell + side.CellOffset;
                if (! isValidCell(face) || ! isValidCell(neighbor))
                    continue;

                int const faceIdx     = index2GridOffset(face);
                int const neighborIdx = index2GridOffset(neighbor);
                float const open      = m_s[neighborIdx] > 0.0f ? 1.0f : 0.0f;
                if (open <= 0.0f)
                    continue;

                divergence += double(side.Sign) * double(m_vel[faceIdx][side.Direction]);

                int const neighborRow = cellToRow[neighborIdx];
                if (neighborRow >= 0) {
                    diagonal += 1.0;
                    if (! pinned[neighborRow])
                        triplets.emplace_back(row, neighborRow, -1.0);
                } else if (m_type[neighborIdx] == EMPTY_CELL) {
                    diagonal += 1.0;
                } else {
                    diagonal += 1.0;
                }
            }

            if (compensateDrift && m_particleRestDensity > 0.0f) {
                float const densityError = m_particleDensity[idx] - m_particleRestDensity;
                if (densityError > 0.0f)
                    divergence -= compensateDriftWeight * densityError;
            }

            if (diagonal <= 1e-8) {
                triplets.emplace_back(row, row, 1.0);
                rhs[row] = 0.0;
            } else {
                triplets.emplace_back(row, row, diagonal);
                rhs[row] = -divergence;
            }
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
        if (solver.info() != Eigen::Success || ! pressure.allFinite())
            pressure.setZero();

        for (int row = 0; row < matrixSize; ++row)
            m_p[rowToCell[row]] = float(pressure[row]);

        for (int row = 0; row < matrixSize; ++row) {
            if (pinned[row])
                continue;

            int const        idx  = rowToCell[row];
            glm::ivec3 const cell = decode(idx);
            float const      p    = m_p[idx];

            for (Neighbor const & side : Neighbors) {
                glm::ivec3 const face     = cell + side.FaceOffset;
                glm::ivec3 const neighbor = cell + side.CellOffset;
                if (! isValidCell(face) || ! isValidCell(neighbor))
                    continue;

                int const neighborIdx = index2GridOffset(neighbor);
                if (m_s[neighborIdx] <= 0.0f)
                    continue;

                int const faceIdx = index2GridOffset(face);
                if (cellToRow[neighborIdx] >= 0 || m_type[neighborIdx] == EMPTY_CELL)
                    m_vel[faceIdx][side.Direction] += side.Sign * p;
            }
        }

        for (glm::vec3 & v : m_vel) {
            if (! finite(v)) {
                v = glm::vec3(0.0f);
                continue;
            }
            float const speed = glm::length(v);
            if (speed > 8.0f)
                v *= 8.0f / speed;
        }

        if (! m_body)
            return;

        float const pressureForceScale = 8.0f;
        for (int row = 0; row < matrixSize; ++row) {
            int const        idx  = rowToCell[row];
            glm::ivec3 const cell = decode(idx);
            glm::vec3 const  pos  = cellCenter(cell);
            float const      dist = m_body->GetSDF(pos);
            float const      press = m_p[idx];

            if (press <= 0.0f || dist >= m_h)
                continue;

            glm::vec3 const normal = m_body->GetSDFNormal(pos, 0.25f * m_h);
            glm::vec3 const force =
                -normal * press * m_h * m_h * pressureForceScale;
            if (! finite(force))
                continue;

            m_feedbackForce += force;
            m_feedbackTorque += glm::cross(pos - m_body->position, force);
        }

        if (! finite(m_feedbackForce))
            m_feedbackForce = glm::vec3(0.0f);
        if (! finite(m_feedbackTorque))
            m_feedbackTorque = glm::vec3(0.0f);
    }
}
