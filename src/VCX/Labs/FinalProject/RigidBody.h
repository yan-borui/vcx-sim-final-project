#pragma once

#ifndef GLM_ENABLE_EXPERIMENTAL
#    define GLM_ENABLE_EXPERIMENTAL
#endif

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

namespace VCX::Labs::Final {
    struct RigidBody {
        enum class ShapeType {
            Box    = 0,
            Sphere = 1,
            Bunny  = 2
        };

        ShapeType shape = ShapeType::Box;

        glm::vec3 position        = { 0, 0, 0 };
        glm::quat orientation     = { 1, 0, 0, 0 };
        glm::vec3 velocity        = { 0, 0, 0 };
        glm::vec3 angularVelocity = { 0, 0, 0 };
        float     mass            = 1.0f;
        glm::mat3 inertiaLocalInv = glm::mat3(1.0f);
        glm::vec3 dim             = { 0.3f, 0.3f, 0.3f };
        glm::vec3 color           = { 0.8f, 0.2f, 0.2f };
        float     restitution     = 0.5f;
        bool      isStatic        = false;

        static float SdSphere(glm::vec3 const & p, float r) {
            return glm::length(p) - r;
        }

        static float SdEllipsoid(glm::vec3 const & p, glm::vec3 const & r) {
            glm::vec3 q  = p / r;
            float     k0 = glm::length(q);
            glm::vec3 rr = r * r;
            float     k1 = glm::length(p / rr);
            return k0 * (k0 - 1.0f) / std::max(k1, 1e-6f);
        }

        static float SdCapsule(glm::vec3 const & p, glm::vec3 const & a, glm::vec3 const & b, float r) {
            glm::vec3 pa = p - a;
            glm::vec3 ba = b - a;
            float     h  = glm::clamp(glm::dot(pa, ba) / std::max(glm::dot(ba, ba), 1e-6f), 0.0f, 1.0f);
            return glm::length(pa - ba * h) - r;
        }

        static float SdBox(glm::vec3 const & p, glm::vec3 const & halfDim) {
            glm::vec3 q           = glm::abs(p) - halfDim;
            float     outsideDist = glm::length(glm::max(q, glm::vec3(0.0f)));
            float     insideDist  = std::min(std::max(q.x, std::max(q.y, q.z)), 0.0f);
            return outsideDist + insideDist;
        }

        static float BunnyLikeSDF(glm::vec3 const & pUnit) {
            // pUnit is roughly in local normalized coordinates around [-0.5, 0.5]
            // Body
            float body = SdEllipsoid(pUnit - glm::vec3(0.0f, -0.02f, 0.0f), glm::vec3(0.22f, 0.16f, 0.18f));

            // Head
            float head = SdSphere(pUnit - glm::vec3(0.0f, 0.16f, 0.10f), 0.10f);

            // Ears
            float earL = SdCapsule(
                pUnit,
                glm::vec3(-0.05f, 0.22f, 0.08f),
                glm::vec3(-0.10f, 0.46f, 0.08f),
                0.035f);

            float earR = SdCapsule(
                pUnit,
                glm::vec3(0.05f, 0.22f, 0.08f),
                glm::vec3(0.10f, 0.48f, 0.08f),
                0.035f);

            // Tail
            float tail = SdSphere(pUnit - glm::vec3(0.0f, 0.02f, -0.18f), 0.06f);

            // Feet
            float footL = SdEllipsoid(pUnit - glm::vec3(-0.08f, -0.16f, 0.08f), glm::vec3(0.07f, 0.04f, 0.10f));
            float footR = SdEllipsoid(pUnit - glm::vec3(0.08f, -0.16f, 0.08f), glm::vec3(0.07f, 0.04f, 0.10f));

            float d = body;
            d       = std::min(d, head);
            d       = std::min(d, earL);
            d       = std::min(d, earR);
            d       = std::min(d, tail);
            d       = std::min(d, footL);
            d       = std::min(d, footR);

            return d;
        }

        void ComputeInertia() {
            if (isStatic || mass <= 0.0f) {
                inertiaLocalInv = glm::mat3(0.0f);
                return;
            }

            if (shape == ShapeType::Sphere) {
                float const r    = 0.5f * dim.x;
                float const I    = 0.4f * mass * r * r;
                float const invI = I > 1e-6f ? 1.0f / I : 0.0f;
                inertiaLocalInv  = glm::mat3(invI);
                return;
            }

            float x2        = dim.x * dim.x;
            float y2        = dim.y * dim.y;
            float z2        = dim.z * dim.z;
            float ix        = (1.0f / 12.0f) * mass * (y2 + z2);
            float iy        = (1.0f / 12.0f) * mass * (x2 + z2);
            float iz        = (1.0f / 12.0f) * mass * (x2 + y2);
            ix              = std::max(ix, 1e-6f);
            iy              = std::max(iy, 1e-6f);
            iz              = std::max(iz, 1e-6f);
            inertiaLocalInv = glm::inverse(glm::mat3(ix, 0, 0, 0, iy, 0, 0, 0, iz));
        }

        glm::mat3 GetInertiaWorldInv() const {
            if (isStatic || mass <= 0.0f)
                return glm::mat3(0.0f);
            glm::mat3 R = glm::mat3_cast(orientation);
            return R * inertiaLocalInv * glm::transpose(R);
        }

        float GetSDF(glm::vec3 const & worldP) const {
            glm::vec3 const relativeP = worldP - position;
            glm::mat3 const invR      = glm::transpose(glm::mat3_cast(orientation));
            glm::vec3 const localP    = invR * relativeP;

            if (shape == ShapeType::Sphere) {
                float const r = 0.5f * dim.x;
                return glm::length(localP) - r;
            }

            if (shape == ShapeType::Bunny) {
                glm::vec3 pUnit = localP / glm::max(dim, glm::vec3(1e-6f));
                return BunnyLikeSDF(pUnit);
            }

            return SdBox(localP, dim * 0.5f);
        }

        void Reset(
            glm::vec3 const & pos,
            glm::vec3 const & vel,
            glm::vec3 const & size,
            float             m,
            glm::vec3 const & col,
            ShapeType         newShape = ShapeType::Box) {
            position        = pos;
            velocity        = vel;
            angularVelocity = { 0, 0, 0 };
            orientation     = { 1, 0, 0, 0 };
            dim             = size;
            mass            = m;
            color           = col;
            shape           = newShape;
            isStatic        = false;
            ComputeInertia();
        }

        glm::vec3 GetVelocityAtPoint(glm::vec3 const & r) const {
            return velocity + glm::cross(angularVelocity, r);
        }

        float BoundingRadius() const {
            if (shape == ShapeType::Sphere)
                return 0.5f * dim.x;
            if (shape == ShapeType::Bunny)
                return 0.5f * std::max(dim.x, std::max(dim.y * 1.3f, dim.z));
            return 0.5f * std::max(dim.x, std::max(dim.y, dim.z));
        }
    };
} // namespace VCX::Labs::Final
