#include "Labs/FinalProject/CGSimulator.h"
#include <Eigen/Sparse>

namespace VCX::Labs::Final {

    void CGSimulator::solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift) {
        // 1. 找出所有流体格子
        std::vector<int> cellToMatId(m_iNumCells, -1);
        std::vector<int> matIdToCell;
        matIdToCell.reserve(m_iNumCells / 4);

        for (int i = 0; i < m_iNumCells; i++) {
            if (m_type[i] == FLUID_CELL) {
                cellToMatId[i] = (int)matIdToCell.size();
                matIdToCell.push_back(i);
            }
        }

        int n = (int)matIdToCell.size();
        if (n == 0) return;

        // 2. 构建线性系统 Ap = b
        Eigen::SparseMatrix<double> A(n, n);
        std::vector<Eigen::Triplet<double>> triplets;
        Eigen::VectorXd b(n);
        triplets.reserve(n * 7);

        for (int id = 0; id < n; id++) {
            int idx = matIdToCell[id];
            int i = idx % m_iCellX;
            int j = (idx / m_iCellX) % m_iCellY;
            int k = idx / (m_iCellX * m_iCellY);

            // 计算散度作为 b
            float div = (m_vel[index2GridOffset({i + 1, j, k})].x - m_vel[idx].x) +
                        (m_vel[index2GridOffset({i, j + 1, k})].y - m_vel[idx].y) +
                        (m_vel[index2GridOffset({i, j, k + 1})].z - m_vel[idx].z);
            
            // IDP 体积补偿
            if (compensateDrift && m_particleRestDensity > 0.0f) {
                float densityError = m_particleDensity[idx] - m_particleRestDensity;
                if (densityError > 0.0f) div -= compensateDriftWeight * densityError;
            }
            b[id] = -div; 

            // 检查 6 个邻居填充矩阵 A
            glm::ivec3 offsets[6] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
            double diag = 0;
            for (auto const& off : offsets) {
                glm::ivec3 neighborPos = glm::ivec3(i,j,k) + off;
                if (neighborPos.x < 0 || neighborPos.x >= m_iCellX || 
                    neighborPos.y < 0 || neighborPos.y >= m_iCellY || 
                    neighborPos.z < 0 || neighborPos.z >= m_iCellZ) continue;

                int nIdx = index2GridOffset(neighborPos);
                if (m_type[nIdx] == FLUID_CELL) {
                    diag += 1.0;
                    triplets.push_back(Eigen::Triplet<double>(id, cellToMatId[nIdx], -1.0));
                } else if (m_type[nIdx] == EMPTY_CELL) {
                    diag += 1.0; // 空气边界 Dirichlet: p=0
                }
            }
            triplets.push_back(Eigen::Triplet<double>(id, id, diag));
        }
        A.setFromTriplets(triplets.begin(), triplets.end());

        // 3. 调用 Eigen CG 求解器
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
        solver.setTolerance(1e-4);
        solver.compute(A);
        Eigen::VectorXd p_sol = solver.solve(b);

        // 4. 将结果写回并更新速度
        std::fill(m_p.begin(), m_p.end(), 0.0f);
        for (int id = 0; id < n; id++) m_p[matIdToCell[id]] = (float)p_sol[id];

        if (m_body) {
            m_feedbackForce  = glm::vec3(0.0f);
            m_feedbackTorque = glm::vec3(0.0f);
            for (int id = 0; id < n; id++) {
                int idx = matIdToCell[id];
                float press = m_p[idx];
                if (press <= 0.0f) continue;

                int i = idx % m_iCellX;
                int j = (idx / m_iCellX) % m_iCellY;
                int k = idx / (m_iCellX * m_iCellY);
                glm::vec3 cellPos = (glm::vec3(i, j, k) + 0.5f) * m_h - 0.5f;

                float dist = m_body->GetSDF(cellPos);
                if (dist < m_h) {
                    float eps = 0.001f;
                    glm::vec3 normal = glm::normalize(glm::vec3(
                        m_body->GetSDF(cellPos + glm::vec3(eps, 0, 0)) - m_body->GetSDF(cellPos - glm::vec3(eps, 0, 0)),
                        m_body->GetSDF(cellPos + glm::vec3(0, eps, 0)) - m_body->GetSDF(cellPos - glm::vec3(0, eps, 0)),
                        m_body->GetSDF(cellPos + glm::vec3(0, 0, eps)) - m_body->GetSDF(cellPos - glm::vec3(0, 0, eps))
                    ));

                    float area = m_h * m_h;
                    glm::vec3 force = -normal * press * area * 45.0f; // 增强反馈力
                    m_feedbackForce += force;

                    glm::vec3 r = cellPos - m_body->position;
                    m_feedbackTorque += glm::cross(r, force);
                }
            }
        }

        for (int id = 0; id < n; id++) {
            int idx = matIdToCell[id];
            int i = idx % m_iCellX; int j = (idx / m_iCellX) % m_iCellY; int k = idx / (m_iCellX * m_iCellY);
            float p = m_p[idx];

            if (m_s[index2GridOffset({i-1, j, k})] > 0) m_vel[idx].x -= p;
            if (m_s[index2GridOffset({i+1, j, k})] > 0) m_vel[index2GridOffset({i+1, j, k})].x += p;
            if (m_s[index2GridOffset({i, j-1, k})] > 0) m_vel[idx].y -= p;
            if (m_s[index2GridOffset({i, j+1, k})] > 0) m_vel[index2GridOffset({i, j+1, k})].y += p;
            if (m_s[index2GridOffset({i, j, k-1})] > 0) m_vel[idx].z -= p;
            if (m_s[index2GridOffset({i, j, k+1})] > 0) m_vel[index2GridOffset({i, j, k+1})].z += p;
        }
    }
}