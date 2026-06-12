#include "APICSimulator.h"

namespace VCX::Labs::Final {
    void APICSimulator::transferVelocities(bool toGrid, float flipRatio) {
        if (toGrid) {
            rebuildSolidCellsFromBody();

            std::fill(m_vel.begin(), m_vel.end(), glm::vec3(0.0f));
            for (int i = 0; i < 3; i++) {
                std::fill(m_near_num[i].begin(), m_near_num[i].end(), 0.0f);
            }

            for (int i = 0; i < m_iNumCells; i++) {
                m_type[i] = (m_s[i] > 0.0f) ? EMPTY_CELL : SOLID_CELL;
            }
        }

        for (int p = 0; p < m_iNumSpheres; p++) {
            glm::vec3 xp = m_particlePos[p];
            glm::mat3 cp = m_particleC[p];
            glm::vec3 vp = m_particleVel[p];

            glm::vec3 posRelGrid = xp + glm::vec3(0.5f);
            glm::ivec3 cellIndex = glm::ivec3(posRelGrid * m_fInvSpacing);
            int centerIdx = index2GridOffset(cellIndex);
            if (centerIdx >= 0 && centerIdx < m_iNumCells && m_type[centerIdx] == EMPTY_CELL) m_type[centerIdx] = FLUID_CELL;   

            // 瀵筙锛孻锛孼涓変釜鏂瑰悜鍒嗗埆澶勭悊
            for (int dir = 0; dir < 3; dir++) {
                glm::vec3 offset = glm::vec3(0.5f);
                offset[dir] = 0.0f; // 鍙湪褰撳墠鏂瑰悜涓婂亸绉?
                glm::vec3 f_idx = (xp + glm::vec3(0.5f) - offset * m_h) * m_fInvSpacing;
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
                    for (int n = 0; n < 8; n++) {
                        int gIdx = index2GridOffset(neighbor[n]);
                        glm::vec3 xi_xp = (glm::vec3(neighbor[n]) + offset) * m_h - posRelGrid;
                        float velocityContribution = vp[dir] + glm::dot(cp[dir], xi_xp);
                        
                        if (gIdx >= 0 && gIdx < m_iNumCells) {
                            m_vel[gIdx][dir] += w[n] * velocityContribution;
                            m_near_num[dir][gIdx] += w[n];
                        }
                    }
                } else {
                    float newV_dir = 0.0f;
                    glm::vec3 newC_dir(0.0f);

                    for (int n = 0; n < 8; n++) {
                        int gIdx = index2GridOffset(neighbor[n]);
                        glm::vec3 xi_xp = (glm::vec3(neighbor[n]) + offset) * m_h - posRelGrid;
                        if (gIdx >= 0 && gIdx < m_iNumCells) {
                            float v_g = m_vel[gIdx][dir];
                            newV_dir += w[n] * v_g;
                            newC_dir += w[n] * v_g * xi_xp;
                        }
                    }

                    m_particleVel[p][dir] = newV_dir; // 鏇存柊绮掑瓙閫熷害
                    m_particleC[p][dir] = newC_dir * (4.0f * m_fInvSpacing * m_fInvSpacing); // 鏇存柊 APIC C
                }
            }
        }

        if (toGrid) {
            // 褰掍竴鍖栫綉鏍奸€熷害
            for (int i = 0; i < m_iNumCells; i++) {
                for (int dir = 0; dir < 3; dir++) {
                    if (m_near_num[dir][i] > 0.0f) {
                        if (isValidVelocity(i % m_iCellX, (i / m_iCellX) % m_iCellY, i / (m_iCellX * m_iCellY), dir)) {
                            m_vel[i][dir] /= m_near_num[dir][i];
                        } else {
                            m_vel[i][dir] = 0.0f; // 娌℃湁鏈夋晥娴佷綋閭诲眳锛岄€熷害璁句负0
                        }
                    }
                }
            }
            m_pre_vel = m_vel; // 淇濆瓨褰撳墠缃戞牸閫熷害鐢ㄤ簬涓嬩竴娆¤绠梫_delta
        }
    }

    void APICSimulator::setupScene(int res) {
        Simulator::setupScene(res);
        m_particleC.clear();
        m_particleC.resize(m_iNumSpheres, glm::mat3(0.0f));
    }
}
