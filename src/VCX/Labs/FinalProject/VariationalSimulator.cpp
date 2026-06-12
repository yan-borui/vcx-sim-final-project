#include "Labs/FinalProject/VariationalSimulator.h"

namespace VCX::Labs::Final {
    void VariationalSimulator::integrateParticles(float timeStep) {
        for (int i = 0; i < m_iNumSpheres; ++i) {
            m_particleVel[i] += gravity * timeStep;
            m_particlePos[i] += m_particleVel[i] * timeStep;
            m_particlePos[i] = glm::clamp(m_particlePos[i], glm::vec3(-0.48f), glm::vec3(0.48f));
        }
    }

    // 计算网格面上的流体权重 (0.0 表示全固体, 1.0 表示全流体)
    float VariationalSimulator::getFaceWeight(int i, int j, int k, int dir) {
        float sumW = 0.0f;
        // 在 0.25 和 0.75 位置采样
        float offsets[2] = { 0.25f * m_h, 0.75f * m_h };

        glm::vec3 basePos = glm::vec3(i, j, k) * m_h - glm::vec3(0.5f);
    
        for (int a = 0; a < 2; ++a) {
            for (int b = 0; b < 2; ++b) {
                glm::vec3 samplePos = basePos;
                if (dir == 0)      samplePos += glm::vec3(0.0f, offsets[a], offsets[b]); // X-face
                else if (dir == 1) samplePos += glm::vec3(offsets[a], 0.0f, offsets[b]); // Y-face
                else               samplePos += glm::vec3(offsets[a], offsets[b], 0.0f); // Z-face

                float dist = m_body->GetSDF(samplePos);
                // 线性平滑
                sumW += std::clamp(0.5f + dist / m_h, 0.0f, 1.0f);
            }
        }
        return sumW * 0.25f; // 取 4 个点的平均值
    }

    void VariationalSimulator::handleParticleCollisions() {
        float buffer = 0.01f * m_h;
        float minBound = -0.5f + m_h;
        float maxBound =  0.5f - m_h;
        for (int i = 0; i < m_iNumSpheres; ++i) {
            for (int d = 0; d < 3; ++d) {
                if (m_particlePos[i][d] < minBound) { 
                    m_particlePos[i][d] = minBound; 
                    if (m_particleVel[i][d] < 0) m_particleVel[i][d] *= -0.0f;
                }
                if (m_particlePos[i][d] > maxBound) { 
                    m_particlePos[i][d] = maxBound; 
                    if (m_particleVel[i][d] > 0) m_particleVel[i][d] *= -0.0f;
                }
            }

            if (m_body) {
                float dist = m_body->GetSDF(m_particlePos[i]);
                if (dist < 0.0f) {
                    // 计算法线 n
                    float eps = 0.001f;
                    glm::vec3 n = glm::normalize(glm::vec3(
                    m_body->GetSDF(m_particlePos[i] + glm::vec3(eps,0,0)) - m_body->GetSDF(m_particlePos[i] - glm::vec3(eps,0,0)),
                    m_body->GetSDF(m_particlePos[i] + glm::vec3(0,eps,0)) - m_body->GetSDF(m_particlePos[i] - glm::vec3(0,eps,0)),
                    m_body->GetSDF(m_particlePos[i] + glm::vec3(0,0,eps)) - m_body->GetSDF(m_particlePos[i] - glm::vec3(0,0,eps))
                    ));

                    // 1. 位置修正
                    m_particlePos[i] -= dist * n; 

                    // 2. 速度修正 
                    glm::vec3 v_rel = m_particleVel[i] - m_body->GetVelocityAtPoint(m_particlePos[i] - m_body->position);
                    float v_normal = glm::dot(v_rel, n);
        
                    if (v_normal < 0) {
                        // 只消除法向相对速度，不影响切向
                        m_particleVel[i] -= 1.0f * v_normal * n;
                    }
                }
            }
        }
    }

    void VariationalSimulator::pushParticlesApart(int numIters) {
        int gridSizeX = m_iCellX;
        int gridSizeY = m_iCellY;
        int gridSizeZ = m_iCellZ;
        int numCells = m_iNumCells;

        std::fill(m_hashtableindex.begin(), m_hashtableindex.end(), 0);
        std::vector<int> particleCellIndex(m_iNumSpheres);

        for (int i = 0; i < m_iNumSpheres; i++) {
            glm::vec3 pos = m_particlePos[i] + glm::vec3(0.5f);
            int cx = std::clamp((int)(pos.x * m_fInvSpacing), 0, gridSizeX - 1);
            int cy = std::clamp((int)(pos.y * m_fInvSpacing), 0, gridSizeY - 1);
            int cz = std::clamp((int)(pos.z * m_fInvSpacing), 0, gridSizeZ - 1);
            int cellIndex = index2GridOffset({ cx, cy, cz });
            particleCellIndex[i] = cellIndex;
            m_hashtableindex[cellIndex + 1]++;
        }

        for (int i = 1; i < numCells; i++) {
            m_hashtableindex[i] += m_hashtableindex[i - 1];
        }

        std::vector<int> currentCellOffset = m_hashtableindex;
        for (int i = 0; i < m_iNumSpheres; i++) {
            int cellIndex = particleCellIndex[i];
            m_hashtable[currentCellOffset[cellIndex]++] = i;
        }

        float minDist = 3.0f * m_particleRadius;
        for (int iter = 0; iter < numIters; iter++) {
            for (int i = 0; i < m_iNumSpheres; i++) {
                glm::vec3 pos = m_particlePos[i] + glm::vec3(0.5f);
                int cx = (int)(pos.x * m_fInvSpacing);
                int cy = (int)(pos.y * m_fInvSpacing);
                int cz = (int)(pos.z * m_fInvSpacing);

                for (int x = cx - 1; x <= cx + 1; x++) {
                    for (int y = cy - 1; y <= cy + 1; y++) {
                        for (int z = cz - 1; z <= cz + 1; z++) {
                            if (x < 0 || x >= gridSizeX || y < 0 || y >= gridSizeY || z < 0 || z >= gridSizeZ) continue;
                            int cellIndex = index2GridOffset({ x, y, z });
                            for (int k = m_hashtableindex[cellIndex]; k < m_hashtableindex[cellIndex + 1]; k++) {
                                int neighborIdx = m_hashtable[k];
                                if (neighborIdx == i) continue;

                                glm::vec3 diff = m_particlePos[i] - m_particlePos[neighborIdx];
                                float dist = glm::length(diff);
                                if (dist < minDist && dist > 0.0f) {
                                    glm::vec3 correction = (minDist - dist) * (diff / dist) * 0.5f;
                                    m_particlePos[i] += correction;
                                    m_particlePos[neighborIdx] -= correction;
                                }
                            }
                        }
                    }
                }
                float wallMin = -0.5f + m_h;
                float wallMax =  0.5f - m_h;

                for (int d = 0; d < 3; d++) {
                    // 检查最小边界 (地面、左墙、后墙)
                    float dMin = m_particlePos[i][d] - wallMin;
                    if (dMin < m_particleRadius) { // 如果粒子半径触碰到了墙
                        m_particlePos[i][d] += (m_particleRadius - dMin) * 0.5f; // 往回推
                    }
                    // 检查最大边界 (天花板、右墙、前墙)
                    float dMax = wallMax - m_particlePos[i][d];
                    if (dMax < m_particleRadius) {
                        m_particlePos[i][d] -= (m_particleRadius - dMax) * 0.5f; // 往回推
                    }
                }
            }
        }
    }

    void VariationalSimulator::updateParticleDensity() {
        std::fill(m_particleDensity.begin(), m_particleDensity.end(), 0.0f);
        for (int p = 0; p < m_iNumSpheres; p++) {
            glm::vec3 pos = m_particlePos[p] + glm::vec3(0.5f);
            glm::ivec3 cellIndex = glm::ivec3(pos * m_fInvSpacing);

            int idx = index2GridOffset(cellIndex);
            if (idx >= 0 && idx < m_iNumCells) {
                m_particleDensity[idx] += 1.0f; // 每个粒子对所在格子的密度贡献为1
            }
        }
        if (m_particleRestDensity == 0.0f) {
            float sum = 0;
            int count = 0;
            for (int i = 0; i < m_iNumCells; i++) {
                sum += m_particleDensity[i];
                count++;
            }
            m_particleRestDensity = sum / count; // 计算平均密度作为静止状态下的密度
        }
    }

    void VariationalSimulator::transferVelocities(bool toGrid, float flipRatio) {
        if (toGrid) {
            // clear grid velocities and near_num
            std::fill(m_vel.begin(), m_vel.end(), glm::vec3(0.0f));
            for (int i = 0; i < 3; i++) {
                std::fill(m_near_num[i].begin(), m_near_num[i].end(), 0.0f);
            }

            for (int i = 0; i < m_iNumCells; i++) {
                m_type[i] = (m_s[i] > 0.0f) ? EMPTY_CELL : SOLID_CELL;
            }
        }

        // 遍历所有粒子
        for (int p = 0; p < m_iNumSpheres; p++) {
            glm::vec3 pos = m_particlePos[p];
            glm::vec3 posRelGrid = pos + glm::vec3(0.5f);
            glm::ivec3 cellIndex = glm::ivec3(posRelGrid * m_fInvSpacing);

            int centerIdx = index2GridOffset(cellIndex);
            if (centerIdx >= 0 && centerIdx < m_iNumCells && m_type[centerIdx] == EMPTY_CELL) m_type[centerIdx] = FLUID_CELL;

            // 对X，Y，Z三个方向分别处理
            for (int dir = 0; dir < 3; dir++) {
                glm::vec3 offset = glm::vec3(0.5f);
                offset[dir] = 0.0f; // 只在当前方向上偏移

                glm::vec3 f_idx = (pos + glm::vec3(0.5f) - offset * m_h) * m_fInvSpacing;
                glm::ivec3 baseIdx = glm::ivec3(floor(f_idx));
                glm::vec3 delta = f_idx - glm::vec3(baseIdx);
                glm::vec3 deltaC = glm::vec3(1.0f) - delta;

                float w[8];
                w[0] = deltaC.x  * deltaC.y * deltaC.z;
                w[1] = delta.x   * deltaC.y * deltaC.z;
                w[2] = deltaC.x  * delta.y  * deltaC.z;
                w[3] = delta.x   * delta.y  * deltaC.z;
                w[4] = deltaC.x  * deltaC.y * delta.z;
                w[5] = delta.x   * deltaC.y * delta.z;
                w[6] = deltaC.x  * delta.y  * delta.z;
                w[7] = delta.x   * delta.y  * delta.z;

                glm::ivec3 neighbor[8] = {
                    baseIdx + glm::ivec3(0, 0, 0), baseIdx + glm::ivec3(1, 0, 0), 
                    baseIdx + glm::ivec3(0, 1, 0), baseIdx + glm::ivec3(1, 1, 0),
                    baseIdx + glm::ivec3(0, 0, 1), baseIdx + glm::ivec3(1, 0, 1), 
                    baseIdx + glm::ivec3(0, 1, 1), baseIdx + glm::ivec3(1, 1, 1)
                };

                if (toGrid) {
                    // P2G: 粒子到网格
                    float pVel = m_particleVel[p][dir];
                    for (int n = 0; n < 8; n++) {
                        int gIdx = index2GridOffset(neighbor[n]);
                        if (gIdx >= 0 && gIdx < m_iNumCells) {
                            m_vel[gIdx][dir] += w[n] * pVel;
                            m_near_num[dir][gIdx] += w[n];
                        }
                    }
                } else {
                    // G2P: 网格到粒子
                    float v_pic = 0.0f;
                    float v_delta = 0.0f;
                    for (int n = 0; n < 8; n++) {
                        int gIdx = index2GridOffset(neighbor[n]);
                        if (gIdx >= 0 && gIdx < m_iNumCells) {
                            v_pic += w[n] * m_vel[gIdx][dir];
                            v_delta += w[n] * (m_vel[gIdx][dir] - m_pre_vel[gIdx][dir]);
                        }
                    }

                    float v_flip = m_particleVel[p][dir] + v_delta;
                    m_particleVel[p][dir] = flipRatio * v_flip + (1.0f - flipRatio) * v_pic;
                }
            }
        }

        if (toGrid) {
            // 归一化网格速度
            for (int i = 0; i < m_iNumCells; i++) {
                for (int dir = 0; dir < 3; dir++) {
                    if (m_near_num[dir][i] > 0.0001f) {
                        if (isValidVelocity(i % m_iCellX, (i / m_iCellX) % m_iCellY, i / (m_iCellX * m_iCellY), dir)) {
                            m_vel[i][dir] /= m_near_num[dir][i];
                            if (std::abs(m_vel[i][dir]) > 5.0f) {
                                m_vel[i][dir] = (m_vel[i][dir] > 0) ? 5.0f : -5.0f;
                            }
                        } else {
                            m_vel[i][dir] = 0.0f; // 没有有效流体邻居，速度设为0
                        }
                    }
                }
            }
            m_pre_vel = m_vel; // 保存当前网格速度用于下一次计算v_delta
        }
    }

    void VariationalSimulator::solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift) {
        int n = m_iNumCells;
        Eigen::SparseMatrix<float> A(n, n);
        Eigen::VectorXf b(n);
        b.setZero();
        std::vector<Eigen::Triplet<float>> triplets;

        float invH = 1.0f / m_h;
        float scale = 1.0f / (m_h * m_h);

        // --- 1: 预计算刚体逆质量矩阵 ---
        Eigen::Matrix<float, 6, 6> invMs;
        invMs.setZero();
        if (m_body && !m_body->isStatic) {
            invMs.block<3, 3>(0, 0) = Eigen::Matrix3f::Identity() * (1.0f / m_body->mass);
            glm::mat3 invI = m_body->GetInertiaWorldInv();
            for(int r = 0; r < 3; ++r) 
                for(int c = 0; c < 3; ++c) 
                    invMs(3+r, 3+c) = invI[c][r];
        }

        // --- 2: 收集 Jacobian 向量 ---
        std::vector<Eigen::Matrix<float, 6, 1>> J(n, Eigen::Matrix<float, 6, 1>::Zero());
        std::vector<int> boundaryCells;

        for (int i = 1; i < m_iCellX - 1; i++) {
            for (int j = 1; j < m_iCellY - 1; j++) {
                for (int k = 1; k < m_iCellZ - 1; k++) {
                    int idx = index2GridOffset({i, j, k});
                    if (m_type[idx] != FLUID_CELL) continue;

                    // 检查 6 个面
                    int fIdx[6][3] = {{i,j,k}, {i+1,j,k}, {i,j,k}, {i,j+1,k}, {i,j,k}, {i,j+1,k}};
                    int dirs[6] = {0, 0, 1, 1, 2, 2};
                    float signs[6] = {-1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f};

                    for (int m = 0; m < 6; m++) {
                        float w = getFaceWeight(fIdx[m][0], fIdx[m][1], fIdx[m][2], dirs[m]);
                        // w < 1 表示这个面有固体覆盖
                        if (w < 1.0f && m_body) {
                            float area = (1.0f - w) * m_h * m_h;
                            glm::vec3 n_vec(0); n_vec[dirs[m]] = signs[m];
                            glm::vec3 facePos = glm::vec3(fIdx[m][0], fIdx[m][1], fIdx[m][2]) * m_h - 0.5f;
                            glm::vec3 arm = facePos - m_body->position;
                            glm::vec3 torque = glm::cross(arm, n_vec);

                            // Ji 代表压力对刚体产生的广义力
                            J[idx](dirs[m]) += n_vec[dirs[m]] * area;
                            J[idx](3) += torque.x * area;
                            J[idx](4) += torque.y * area;
                            J[idx](5) += torque.z * area;
                        }
                    }
                    if (J[idx].norm() > 1e-9f) boundaryCells.push_back(idx);
                }
            }
        }

        // --- 3: 构建标准流体拉普拉斯矩阵 ---
        for (int i = 1; i < m_iCellX - 1; i++) {
            for (int j = 1; j < m_iCellY - 1; j++) {
                for (int k = 1; k < m_iCellZ - 1; k++) {
                    int idx = index2GridOffset({i, j, k});
                    if (m_type[idx] != FLUID_CELL) continue;

                    float centerDiag = 0.0f;
                    float divergence = 0.0f;
                    float eps_diag = 1e-4f;
                    int neighbors[6][3] = {{i-1,j,k}, {i+1,j,k}, {i,j-1,k}, {i,j+1,k}, {i,j,k-1}, {i,j,k+1}};
                    int dirs[6] = {0, 0, 1, 1, 2, 2};
                    int fIdx[6][3] = {{i,j,k}, {i+1,j,k}, {i,j,k}, {i,j+1,k}, {i,j,k}, {i,j+1,k}};

                    for (int m = 0; m < 6; m++) {
                        float w_raw = getFaceWeight(fIdx[m][0], fIdx[m][1], fIdx[m][2], dirs[m]);
                        float w = (w_raw < 0.05f) ? 0.0f : ((w_raw > 0.95f) ? 1.0f : w_raw);
                        if (w <= 0.0f) continue; 

                        int nIdx = index2GridOffset({neighbors[m][0], neighbors[m][1], neighbors[m][2]});
                        float sgn = (m % 2 == 0) ? 1.0f : -1.0f;

                        // RHS: 散度计算
                        float u_star = m_vel[index2GridOffset({fIdx[m][0], fIdx[m][1], fIdx[m][2]})][dirs[m]];
                        float v_solid_star = 0.0f;
                        if (m_body) {
                            glm::vec3 facePos = glm::vec3(fIdx[m][0], fIdx[m][1], fIdx[m][2]) * m_h - 0.5f;
                            v_solid_star = m_body->GetVelocityAtPoint(facePos - m_body->position)[dirs[m]];
                        }
                        divergence += sgn * (w * u_star + (1.0f - w) * v_solid_star) * invH;

                        // Matrix: 流体部分
                        if (m_type[nIdx] == FLUID_CELL) {
                            triplets.push_back(Eigen::Triplet<float>(idx, nIdx, -w * scale));
                            centerDiag += w * scale;
                        } else if (m_type[nIdx] == EMPTY_CELL) {
                            centerDiag += w * scale; // 自由表面边界
                        }
                    }
                    b(idx) = divergence / dt;
                    triplets.push_back(Eigen::Triplet<float>(idx, idx, centerDiag + eps_diag));
                }
            }
        }

        // --- 4: 耦合项 J * Minv * J^T ---
        if (m_body && !m_body->isStatic) {   
            std::vector<Eigen::Matrix<float, 6, 1>> MinvJ(n);
            for (int idx : boundaryCells) {
                MinvJ[idx] = invMs * J[idx];
            }

            // 填充耦合项
            for (size_t a = 0; a < boundaryCells.size(); ++a) {
                int idx_a = boundaryCells[a];
                for (size_t b = a; b < boundaryCells.size(); ++b) {
                    int idx_b = boundaryCells[b];
                
                    // 计算耦合值: Ji^T * (M^-1 * Jj)
                    float val = J[idx_a].dot(MinvJ[idx_b]);
                
                    if (std::abs(val) > 1e-12f) {
                        triplets.push_back(Eigen::Triplet<float>(idx_a, idx_b, val));
                        if (idx_a != idx_b) {
                            triplets.push_back(Eigen::Triplet<float>(idx_b, idx_a, val)); // 保持矩阵对称
                        }
                    }
                }
            }
        }

        // --- 5: 求解器设置 ---
        A.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower|Eigen::Upper> solver;
        solver.setMaxIterations(200); 
        solver.setTolerance(1e-4);
        solver.compute(A);
        Eigen::VectorXf p_sol = solver.solve(b);

        // 处理 NaN 和 Inf，确保压力值在合理范围内
        for (int i = 0; i < m_iNumCells; i++) {
            float p = (i < p_sol.size()) ? p_sol[i] : 0.0f;
            if (!std::isfinite(p)) p = 0.0f; 
            m_p[i] = std::clamp(p, 0.0f, 50.0f);
        }

        // 更新流体网格速度
        for (int i = 1; i < m_iCellX - 1; i++) {
            for (int j = 1; j < m_iCellY - 1; j++) {
                for (int k = 1; k < m_iCellZ - 1; k++) {
                    int idx = index2GridOffset({i, j, k});
                    int neighbors[3][3] = {{i-1,j,k}, {i,j-1,k}, {i,j,k-1}};
                    for (int dir = 0; dir < 3; dir++) {
                        float w = getFaceWeight(i, j, k, dir);
                        if (w > 0.01f) { 
                            int pIdx = index2GridOffset({neighbors[dir][0], neighbors[dir][1], neighbors[dir][2]});
                            
                            // --- 限制速度增量 ---
                            float deltaV = dt * (m_p[idx] - m_p[pIdx]) * invH;
                            float maxV = 0.5f * m_h / dt;
                            deltaV = std::clamp(deltaV, -maxV, maxV);
                            
                            m_vel[idx][dir] -= deltaV;
                        }
                    }
                }
            }
        }

        // 更新刚体速度: V = V* + dt * M^-1 * sum(Ji * pi)
        if (m_body && !m_body->isStatic) {
            Eigen::Matrix<float, 6, 1> totalImpulse = Eigen::Matrix<float, 6, 1>::Zero();
            for (int idx : boundaryCells) {
                totalImpulse += J[idx] * m_p[idx];
            }
        
            Eigen::Matrix<float, 6, 1> deltaV_6 = invMs * totalImpulse;
        
            // --- 线速度更新 ---
            glm::vec3 dv = glm::vec3(deltaV_6(0), deltaV_6(1), deltaV_6(2));
            float maxDV = 5.0f * m_h / dt;
            if (glm::length(dv) > maxDV) dv = glm::normalize(dv) * maxDV;
            m_body->velocity += dv;
                m_body->velocity *= 0.995f;

            // --- 角速度更新 ---
            glm::vec3 dw = glm::vec3(deltaV_6(3), deltaV_6(4), deltaV_6(5));
            float maxDW = 2.0f; // 限制单步旋转弧度
            if (glm::length(dw) > maxDW) dw = glm::normalize(dw) * maxDW;
            m_body->angularVelocity = (m_body->angularVelocity + dw) * 0.95f; 
        
            if (!std::isfinite(m_body->velocity.x)) m_body->velocity = glm::vec3(0);
            if (!std::isfinite(m_body->angularVelocity.x)) m_body->angularVelocity = glm::vec3(0);
        }
    }

    void VariationalSimulator::updateParticleColors() {
        const glm::vec3 defaultColor(0.2f, 0.6f, 0.9f);
        const glm::vec3 purple(0.7f, 0.1f, 0.9f);
        const glm::vec3 yellow(0.9f, 0.9f, 0.2f);
        const glm::vec3 red(0.9f, 0.2f, 0.2f);

        for (int i = 0; i < m_iNumSpheres; i++) {
            if (m_colorMode == ColorMode::Default) {
                m_particleColor[i] = defaultColor;
                continue; 
            } else if (m_colorMode == ColorMode::Velocity) { 
                // 速度模式：蓝 -> 紫
                float speed = glm::length(m_particleVel[i]);
                float ratio = std::clamp(speed / 0.8f, 0.0f, 1.0f);
                m_particleColor[i] = glm::mix(defaultColor, purple, ratio);
            } 
            else {
                // 密度与压强模式需要查询网格数据
                glm::vec3 posRel = (m_particlePos[i] + glm::vec3(0.5f)) * m_fInvSpacing;
                int cx = std::clamp((int)posRel.x, 0, m_iCellX - 1);
                int cy = std::clamp((int)posRel.y, 0, m_iCellY - 1);
                int cz = std::clamp((int)posRel.z, 0, m_iCellZ - 1);
                int idx = index2GridOffset({cx, cy, cz});

                if (m_colorMode == ColorMode::Density) {
                    // 密度：蓝 -> 黄
                    float diff = m_particleDensity[idx] - m_particleRestDensity;
                    float ratio = std::clamp(diff / 2.0f + 0.1f, 0.0f, 1.0f);
                    m_particleColor[i] = glm::mix(defaultColor, yellow, ratio);
                } else if (m_colorMode == ColorMode::Pressure) {
                    // 压强：蓝 -> 红
                    float pressure = m_p[idx];
                    float ratio = std::clamp(pressure / 1.2f, 0.0f, 1.0f);
                    m_particleColor[i] = glm::mix(defaultColor, red, ratio);
                }
            }
        }
    }
}
