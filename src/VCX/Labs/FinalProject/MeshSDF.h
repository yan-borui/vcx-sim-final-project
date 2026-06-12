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
        glm::vec3 const & MinBounds() const { return _minBounds; }
        glm::vec3 const & MaxBounds() const { return _maxBounds; }

    private:
        std::vector<MeshTriangle> _triangles;
        glm::vec3                 _minBounds { 0.0f };
        glm::vec3                 _maxBounds { 0.0f };
        glm::vec3                 _gridMin { 0.0f };
        glm::vec3                 _gridMax { 0.0f };
        glm::vec3                 _gridDx { 1.0f };
        int                       _gridResolution { 0 };
        std::vector<float>        _distanceGrid;

        static float PointTriangleDistance(glm::vec3 const & p, MeshTriangle const & tri);
        static bool  RayIntersectsTriangle(
             glm::vec3 const &    orig,
             glm::vec3 const &    dir,
             MeshTriangle const & tri);

        int   GridOffset(int x, int y, int z) const;
        float ExactSignedDistance(glm::vec3 const & localPoint) const;
        float OutsideBoundsDistance(glm::vec3 const & localPoint) const;
        void  BuildDistanceGrid();
        bool IsInside(glm::vec3 const & p) const;
    };

} // namespace VCX::Labs::Final
