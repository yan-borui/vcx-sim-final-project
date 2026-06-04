#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/matrix_transform.hpp>

namespace VCX::Labs::Final {
    struct RigidBody {
        glm::vec3 position = {0, 0, 0};
        glm::quat orientation = {1, 0, 0, 0};
        glm::vec3 velocity = {0, 0, 0};
        glm::vec3 angularVelocity = {0, 0, 0};
        float mass = 1.0f;
        glm::mat3 inertiaLocalInv;
        glm::vec3 dim = {0.3f, 0.3f, 0.3f};
        glm::vec3 color = {0.8f, 0.2f, 0.2f};
        float restitution = 0.5f;
        bool isStatic = false;

        void ComputeInertia() {
            if (isStatic || mass <= 0) {
                inertiaLocalInv = glm::mat3(0.0f);
                return;
            }
            float x2 = dim.x * dim.x; float y2 = dim.y * dim.y; float z2 = dim.z * dim.z;
            float ix = (1.0f / 12.0f) * mass * (y2 + z2);
            float iy = (1.0f / 12.0f) * mass * (x2 + z2);
            float iz = (1.0f / 12.0f) * mass * (x2 + y2);
            if (ix < 1e-6f) ix = 1e-6f; if (iy < 1e-6f) iy = 1e-6f; if (iz < 1e-6f) iz = 1e-6f;
            inertiaLocalInv = glm::inverse(glm::mat3(ix, 0, 0, 0, iy, 0, 0, 0, iz));
        }

        glm::mat3 GetInertiaWorldInv() const {
            if (isStatic || mass <= 0) return glm::mat3(0.0f);
            glm::mat3 R = glm::mat3_cast(orientation);
            return R * inertiaLocalInv * glm::transpose(R);
        }

        float GetSDF(const glm::vec3& worldP) const {
            // 1. 将点变换到刚体的局部坐标系（考虑旋转和位移）
            glm::vec3 relativeP = worldP - position;
            glm::mat3 invR = glm::transpose(glm::mat3_cast(orientation));
            glm::vec3 localP = invR * relativeP;

            // 2. 在局部坐标系下计算到 AABB 的距离
            glm::vec3 q = glm::abs(localP) - dim * 0.5f;
            float outsideDist = glm::length(glm::max(q, 0.0f));
            float insideDist = std::min(std::max(q.x, std::max(q.y, q.z)), 0.0f);
            
            return outsideDist + insideDist;
        }

        void Reset(const glm::vec3& pos, const glm::vec3& vel, const glm::vec3& size, float m, const glm::vec3& col) {
            position = pos; 
            velocity = vel; 
            angularVelocity = {0, 0, 0}; 
            orientation = {1, 0, 0, 0};
            dim = size; 
            mass = m; 
            color = col; 
            isStatic = false;
            ComputeInertia();
        }

        glm::vec3 GetVelocityAtPoint(const glm::vec3& r) const {
            return velocity + glm::cross(angularVelocity, r);
        }
    };
}