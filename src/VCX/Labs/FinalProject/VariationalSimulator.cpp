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
        // 采样面中心的坐标
        glm::vec3 facePos = glm::vec3(i, j, k) * m_h - glm::vec3(0.5f);
        if (dir == 0) facePos += glm::vec3(0.0f, 0.5f, 0.5f) * m_h;         // X-face
        else if (dir == 1) facePos += glm::vec3(0.5f, 0.0f, 0.5f) * m_h;    // Y-face
        else facePos += glm::vec3(0.5f, 0.5f, 0.0f) * m_h;                  // Z-face

        float dist = m_body->GetSDF(facePos);
    
        return std::clamp(0.5f + dist / m_h, 0.0f, 1.0f);
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
                        m_particleVel[i] -= v_normal * n;
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

        // --- 1. 构建全局稀疏矩阵与 RHS (b) ---
        for (int i = 1; i < m_iCellX - 1; i++) {
            for (int j = 1; j < m_iCellY - 1; j++) {
                for (int k = 1; k < m_iCellZ - 1; k++) {
                    int idx = index2GridOffset({i, j, k});
                    if (m_type[idx] != FLUID_CELL) continue;

                    float centerDiag = 0.0f;
                    float divergence = 0.0f; 
                
                    int neighbors[6][3] = {{i-1,j,k}, {i+1,j,k}, {i,j-1,k}, {i,j+1,k}, {i,j,k-1}, {i,j,k+1}};
                    int dirs[6] = {0, 0, 1, 1, 2, 2};
                    int faceIdx[6][3] = {{i,j,k}, {i+1,j,k}, {i,j,k}, {i,j+1,k}, {i,j,k}, {i,j,k+1}};

                    for (int m = 0; m < 6; m++) {
                        float w = getFaceWeight(faceIdx[m][0], faceIdx[m][1], faceIdx[m][2], dirs[m]);

                        if (w < 0.05f) w = 0.0f; 
                        else if (w > 0.99f) w = 1.0f;

                        int nIdx = index2GridOffset({neighbors[m][0], neighbors[m][1], neighbors[m][2]});
                        int fGrid = index2GridOffset({faceIdx[m][0], faceIdx[m][1], faceIdx[m][2]});
                        float sgn = (m % 2 == 0) ? 1.0f : -1.0f;

                        // 计算面上的流量贡献
                        float u_star = m_vel[fGrid][dirs[m]];
                        float v_s = 0.0f;
                    
                        if (m_type[nIdx] == SOLID_CELL) {
                            // 如果邻居是固体，获取固体的边界速度
                            glm::vec3 fPos = (glm::vec3(faceIdx[m][0], faceIdx[m][1], faceIdx[m][2])) * m_h - 0.5f;
                            if (m_body && m_body->GetSDF(fPos) < m_h) {
                                v_s = m_body->GetVelocityAtPoint(fPos - m_body->position)[dirs[m]];
                            }
                        }

                        // 散度计算：w*u + (1-w)*v_s
                        divergence += sgn * (w * u_star + (1.0f - w) * v_s) * invH;

                        if (m_type[nIdx] == FLUID_CELL) {
                            triplets.push_back(Eigen::Triplet<float>(idx, nIdx, -w * scale));
                            centerDiag += w * scale;
                        } else if (m_type[nIdx] == EMPTY_CELL) {
                            centerDiag += w * scale;
                        }
                    }

                    b(idx) = divergence / dt;

                    if (compensateDrift && m_particleRestDensity > 0.001f) {
                        float targetDensity = m_particleRestDensity * 0.8f; 
                        float dErr = m_particleDensity[idx] - targetDensity;
    
                        if (dErr > 0) {
                            b(idx) += 0.1f * (dErr / dt); 
                        }
                    }
                    triplets.push_back(Eigen::Triplet<float>(idx, idx, centerDiag));
                }
            }
        }

        // --- 2. 求解 ---
        for (int i = 0; i < b.size(); ++i) {
            if (std::isnan(b(i))) b(i) = 0.0f;
            b(i) = std::clamp(b(i), -20.0f, 20.0f); 
        }
    

        A.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower|Eigen::Upper> solver;
        solver.setMaxIterations(150); 
        solver.setTolerance(0.01f); 

        solver.compute(A);
        Eigen::VectorXf p_sol = solver.solve(b);

        // 限制压力为非负，壁面分离逻辑
        for (int i = 0; i < p_sol.size(); i++) {
        if (p_sol[i] < 0.0f || std::isnan(p_sol[i])) p_sol[i] = 0.0f;
    
        // 限制压强上限
        if (p_sol[i] > 50.0f) p_sol[i] = 50.0f; 
            m_p[i] = p_sol[i];
        }

        // --- 3. 同步更新速度与反馈冲量 ---
        m_feedbackForce = glm::vec3(0.0f);
        m_feedbackTorque = glm::vec3(0.0f);

        for (int i = 1; i < m_iCellX - 1; i++) {
            for (int j = 1; j < m_iCellY - 1; j++) {
                for (int k = 1; k < m_iCellZ - 1; k++) {
                    int idx = index2GridOffset({i, j, k});

                    float p_val = m_p[idx];
                    if (p_val <= 0.0f) continue;

                    int neighbors[6][3] = {{i-1,j,k}, {i+1,j,k}, {i,j-1,k}, {i,j+1,k}, {i,j,k-1}, {i,j,k+1}};
                    int dirs[6] = {0, 0, 1, 1, 2, 2};
                    // 面坐标对齐：负向面(0,2,4)就是(i,j,k)，正向面(1,3,5)是(i+1,j+1,k+1)
                    int faces[6][3] = {{i,j,k}, {i+1,j,k}, {i,j,k}, {i,j+1,k}, {i,j,k}, {i,j,k+1}};
                    float signs[6] = {-1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f}; // 面的法线方向

                    for (int m = 0; m < 6; m++) {
                        float w = getFaceWeight(faces[m][0], faces[m][1], faces[m][2], dirs[m]);
                
                        // 如果面上有固体成分 (w < 1.0)
                        if (w < 0.999f && m_body) {
                            // 压强产生的力：F = P * Area * (1-w)
                            // 冲量 J = F * dt
                            float solidAreaFrac = (1.0f - w);
                            float impulseMag = p_val * solidAreaFrac * (m_h * m_h) * dt;
                    
                            glm::vec3 n(0); n[dirs[m]] = signs[m];
                            glm::vec3 jVec = n * impulseMag;
                    
                            // 采样点：面的中心
                            glm::vec3 facePos = glm::vec3(faces[m][0], faces[m][1], faces[m][2]) * m_h - 0.5f;
                            glm::vec3 offset(0.5f * m_h); offset[dirs[m]] = 0.0f;
                            glm::vec3 actualPos = facePos + offset;

                            m_feedbackForce  += jVec;
                            m_feedbackTorque += glm::cross(actualPos - m_body->position, jVec);
                        }
                    }
                }
            }
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
