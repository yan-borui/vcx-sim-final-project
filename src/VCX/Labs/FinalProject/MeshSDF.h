#pragma once

#include <string>
#include <vector>

#include <glm/glm.hpp>

namespace VCX::Labs::Final {

    struct MeshTriangle {
        glm::vec3 a;
        glm::vec3 b;
        glm::vec3 c;
    };

    class MeshSDF {
    public:
        bool  LoadOBJ(std::string const & filePath);
        float SignedDistance(glm::vec3 const & localPoint) const;
        bool  IsValid() const { return ! _triangles.empty(); }

    private:
        std::vector<MeshTriangle> _triangles;

        static float PointTriangleDistance(glm::vec3 const & p, MeshTriangle const & tri);
        static bool  RayIntersectsTriangle(
             glm::vec3 const &    orig,
             glm::vec3 const &    dir,
             MeshTriangle const & tri);

        bool IsInside(glm::vec3 const & p) const;
    };

} // namespace VCX::Labs::Final