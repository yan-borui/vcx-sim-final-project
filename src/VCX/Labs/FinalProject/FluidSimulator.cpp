#include "Labs/FinalProject/FluidSimulator.h"

namespace VCX::Labs::Final {
    std::vector<float> Simulator::buildParticleLevelSet() const {
        float const farDistance         = 3.0f * m_h;
        float const reconstructionRadius = 0.75f * m_h;
        std::vector<float> levelSet(m_iNumCells, farDistance);

        auto offset = [&](glm::ivec3 const cell) {
            return cell.x + m_iCellX * (cell.y + m_iCellY * cell.z);
        };

        for (glm::vec3 const particle : m_particlePos) {
            glm::ivec3 const base = glm::clamp(
                glm::ivec3((particle + glm::vec3(0.5f)) * m_fInvSpacing),
                glm::ivec3(0),
                glm::ivec3(m_iCellX - 1, m_iCellY - 1, m_iCellZ - 1));
            for (int z = -2; z <= 2; ++z) {
                for (int y = -2; y <= 2; ++y) {
                    for (int x = -2; x <= 2; ++x) {
                        glm::ivec3 const cell = base + glm::ivec3(x, y, z);
                        if (cell.x < 0 || cell.x >= m_iCellX
                            || cell.y < 0 || cell.y >= m_iCellY
                            || cell.z < 0 || cell.z >= m_iCellZ)
                            continue;

                        glm::vec3 const center =
                            (glm::vec3(cell) + glm::vec3(0.5f)) * m_h
                            - glm::vec3(0.5f);
                        int const idx = offset(cell);
                        levelSet[idx] = std::min(
                            levelSet[idx],
                            glm::length(center - particle) - reconstructionRadius);
                    }
                }
            }
        }

        for (int idx = 0; idx < m_iNumCells; ++idx) {
            if (m_type[idx] == FLUID_CELL && levelSet[idx] >= 0.0f)
                levelSet[idx] = -0.5f * m_h;
            else if (m_type[idx] == EMPTY_CELL && levelSet[idx] <= 0.0f)
                levelSet[idx] = 0.5f * m_h;
        }
        return levelSet;
    }

    float Simulator::ghostFluidPressureScale(
        std::vector<float> const & levelSet,
        int                        fluidCell,
        int                        airCell) const {
        if (fluidCell < 0 || fluidCell >= int(levelSet.size())
            || airCell < 0 || airCell >= int(levelSet.size()))
            return 2.0f;

        float const minimumDistance = 1e-4f * m_h;
        float const fluidDistance =
            std::min(levelSet[fluidCell], -minimumDistance);
        float const airDistance =
            std::max(levelSet[airCell], minimumDistance);
        float const theta = std::clamp(
            -fluidDistance / (airDistance - fluidDistance),
            0.1f,
            1.0f);
        return 1.0f / theta;
    }

    void Simulator::integrateParticles(float timeStep) {
        for (int i = 0; i < m_iNumSpheres; ++i) {
            m_particleVel[i] += gravity * timeStep;
            m_particlePos[i] += m_particleVel[i] * timeStep;
        }
    }

    void Simulator::handleParticleCollisions() {
        float minBound = -0.5f + m_h;
        float maxBound =  0.5f - m_h;
        for (int i = 0; i < m_iNumSpheres; ++i) {
            for (int d = 0; d < 3; ++d) {
                if (m_particlePos[i][d] < minBound) { 
                    m_particlePos[i][d] = minBound; 
                    if (m_particleVel[i][d] < 0) m_particleVel[i][d] *= -0.1f;
                }
                if (m_particlePos[i][d] > maxBound) { 
                    m_particlePos[i][d] = maxBound; 
                    if (m_particleVel[i][d] > 0) m_particleVel[i][d] *= -0.1f;
                }
            }

            if (m_body) {
                float dist = m_body->GetSDF(m_particlePos[i]);
                if (dist < 0.0f) {
                    // 计算表面法线 (通过微扰采样)
                    float eps = 0.001f;
                    glm::vec3 n = glm::normalize(glm::vec3(
                        m_body->GetSDF(m_particlePos[i] + glm::vec3(eps,0,0)) - m_body->GetSDF(m_particlePos[i] - glm::vec3(eps,0,0)),
                        m_body->GetSDF(m_particlePos[i] + glm::vec3(0,eps,0)) - m_body->GetSDF(m_particlePos[i] - glm::vec3(0,eps,0)),
                        m_body->GetSDF(m_particlePos[i] + glm::vec3(0,0,eps)) - m_body->GetSDF(m_particlePos[i] - glm::vec3(0,0,eps))
                    ));
                    m_particlePos[i] -= dist * n; // 推回表面
                    glm::vec3 const bodyVelocity =
                        m_body->GetVelocityAtPoint(m_particlePos[i] - m_body->position);
                    float const inwardSpeed =
                        glm::dot(m_particleVel[i] - bodyVelocity, n);
                    if (inwardSpeed < 0.0f)
                        m_particleVel[i] -= inwardSpeed * n;
                }
            }
        }
    }

    void Simulator::pushParticlesApart(int numIters) {
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

        float minDist = 2.0f * m_particleRadius;
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

    void Simulator::updateParticleDensity() {
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

    void Simulator::transferVelocities(bool toGrid, float flipRatio) {
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
                    if (m_near_num[dir][i] > 0.0f) {
                        if (isValidVelocity(i % m_iCellX, (i / m_iCellX) % m_iCellY, i / (m_iCellX * m_iCellY), dir)) {
                            m_vel[i][dir] /= m_near_num[dir][i];
                        } else {
                            m_vel[i][dir] = 0.0f; // 没有有效流体邻居，速度设为0
                        }
                    }
                }
            }
            m_pre_vel = m_vel; // 保存当前网格速度用于下一次计算v_delta
        }
    }
    
    void Simulator::solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift) {
        std::fill(m_p.begin(), m_p.end(), 0.0f);

        for (int iter = 0; iter < numIters; iter++) {
            for (int i = 1; i < m_iCellX - 1; i++) {
                for (int j = 1; j < m_iCellY - 1; j++) {
                    for (int k = 1; k < m_iCellZ - 1; k++) {
                        int idx = index2GridOffset({ i, j, k });
                        if (m_type[idx] == FLUID_CELL) {
                            float s_left    = m_s[index2GridOffset({ i - 1, j, k })];
                            float s_right   = m_s[index2GridOffset({ i + 1, j, k })];
                            float s_down    = m_s[index2GridOffset({ i, j - 1, k })];
                            float s_up      = m_s[index2GridOffset({ i, j + 1, k })];
                            float s_back    = m_s[index2GridOffset({ i, j, k - 1 })];
                            float s_front   = m_s[index2GridOffset({ i, j, k + 1 })];

                            float s_total = s_left + s_right + s_down + s_up + s_back + s_front;
                            if (s_total > 0.0f) {
                                float divergence = (m_vel[index2GridOffset({ i + 1, j, k })].x - m_vel[idx].x) +
                                                   (m_vel[index2GridOffset({ i, j + 1, k })].y - m_vel[idx].y) +
                                                   (m_vel[index2GridOffset({ i, j, k + 1 })].z - m_vel[idx].z);

                                if (compensateDrift && m_particleRestDensity > 0.0f) {
                                    float densityError = m_particleDensity[idx] - m_particleRestDensity;
                                    if (densityError > 0.0f) {
                                        divergence -= compensateDriftWeight * densityError;
                                    }
                                }

                                float p = (-divergence / s_total) * overRelaxation;
                                m_p[idx] += p;

                                m_vel[idx].x -= s_left * p;
                                m_vel[index2GridOffset({ i + 1, j, k })].x += s_right * p;
                                m_vel[idx].y -= s_down * p;
                                m_vel[index2GridOffset({ i, j + 1, k })].y += s_up * p;
                                m_vel[idx].z -= s_back * p;
                                m_vel[index2GridOffset({ i, j, k + 1 })].z += s_front * p;

                            }
                        }
                    }
                }
            }
        }

        if (m_body) {
            m_feedbackForce  = glm::vec3(0.0f);
            m_feedbackTorque = glm::vec3(0.0f);

            for (int i = 1; i < m_iCellX - 1; i++) {
                for (int j = 1; j < m_iCellY - 1; j++) {
                    for (int k = 1; k < m_iCellZ - 1; k++) {
                        int idx = index2GridOffset({i, j, k});
                        float press = m_p[idx];
                        if (press <= 0.0f) continue;

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
                            glm::vec3 force = -normal * press * area * 45.0f;
                            m_feedbackForce += force;

                            glm::vec3 r = cellPos - m_body->position;
                            m_feedbackTorque += glm::cross(r, force);
                        }
                    }
                }
            }
        }
    }
    
    void Simulator::updateParticleColors() {
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
