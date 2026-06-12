#include "Labs/FinalProject/MeshSDF.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#include <tiny_obj_loader.h>

namespace VCX::Labs::Final {

    namespace {
        float MaxComponent(glm::vec3 const & v) {
            return std::max(v.x, std::max(v.y, v.z));
        }

        glm::vec3 ClosestPointOnTriangle(
            glm::vec3 const & p,
            glm::vec3 const & a,
            glm::vec3 const & b,
            glm::vec3 const & c) {
            glm::vec3 ab = b - a;
            glm::vec3 ac = c - a;
            glm::vec3 ap = p - a;

            float d1 = glm::dot(ab, ap);
            float d2 = glm::dot(ac, ap);
            if (d1 <= 0.0f && d2 <= 0.0f) return a;

            glm::vec3 bp = p - b;
            float     d3 = glm::dot(ab, bp);
            float     d4 = glm::dot(ac, bp);
            if (d3 >= 0.0f && d4 <= d3) return b;

            float vc = d1 * d4 - d3 * d2;
            if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
                float v = d1 / (d1 - d3);
                return a + v * ab;
            }

            glm::vec3 cp = p - c;
            float     d5 = glm::dot(ab, cp);
            float     d6 = glm::dot(ac, cp);
            if (d6 >= 0.0f && d5 <= d6) return c;

            float vb = d5 * d2 - d1 * d6;
            if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
                float w = d2 / (d2 - d6);
                return a + w * ac;
            }

            float va = d3 * d6 - d5 * d4;
            if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
                float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
                return b + w * (c - b);
            }

            float denom = 1.0f / (va + vb + vc);
            float v     = vb * denom;
            float w     = vc * denom;
            return a + ab * v + ac * w;
        }
    } // namespace

    bool MeshSDF::LoadOBJ(std::string const & filePath) {
        std::ifstream ifs(filePath, std::ios::binary);
        if (! ifs.is_open()) {
            std::cerr << "ifstream cannot open file: " << filePath << std::endl;
            return false;
        }

        std::stringstream buffer;
        buffer << ifs.rdbuf();
        std::string objText = buffer.str();
        ifs.close();

        tinyobj::ObjReader       reader;
        tinyobj::ObjReaderConfig config;
        config.triangulate = true;

        if (! reader.ParseFromString(objText, "", config)) {
            if (! reader.Error().empty()) {
                std::cerr << "tinyobj error: " << reader.Error() << std::endl;
            }
            return false;
        }

        if (! reader.Warning().empty()) {
            std::cerr << "tinyobj warning: " << reader.Warning() << std::endl;
        }

        auto const & attrib = reader.GetAttrib();
        auto const & shapes = reader.GetShapes();

        _triangles.clear();
        _minBounds = glm::vec3(std::numeric_limits<float>::max());
        _maxBounds = glm::vec3(std::numeric_limits<float>::lowest());

        for (auto const & shape : shapes) {
            size_t indexOffset = 0;
            for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); ++f) {
                int fv = shape.mesh.num_face_vertices[f];
                if (fv != 3) {
                    indexOffset += size_t(fv);
                    continue;
                }

                glm::vec3 verts[3];
                for (int v = 0; v < 3; ++v) {
                    tinyobj::index_t idx = shape.mesh.indices[indexOffset + size_t(v)];
                    int              vi  = idx.vertex_index;
                    verts[v]             = glm::vec3(
                        attrib.vertices[3 * vi + 0],
                        attrib.vertices[3 * vi + 1],
                        attrib.vertices[3 * vi + 2]);
                    _minBounds = glm::min(_minBounds, verts[v]);
                    _maxBounds = glm::max(_maxBounds, verts[v]);
                }

                _triangles.push_back({ verts[0], verts[1], verts[2] });
                indexOffset += 3;
            }
        }

        if (_triangles.empty())
            return false;

        glm::vec3 const center = 0.5f * (_minBounds + _maxBounds);
        float const     scale  =
            1.0f / std::max(MaxComponent(_maxBounds - _minBounds), 1e-6f);

        _minBounds = glm::vec3(std::numeric_limits<float>::max());
        _maxBounds = glm::vec3(std::numeric_limits<float>::lowest());
        for (MeshTriangle & tri : _triangles) {
            tri.a = (tri.a - center) * scale;
            tri.b = (tri.b - center) * scale;
            tri.c = (tri.c - center) * scale;
            _minBounds = glm::min(_minBounds, glm::min(tri.a, glm::min(tri.b, tri.c)));
            _maxBounds = glm::max(_maxBounds, glm::max(tri.a, glm::max(tri.b, tri.c)));
        }

        BuildDistanceGrid();

        return true;
    }

    float MeshSDF::PointTriangleDistance(glm::vec3 const & p, MeshTriangle const & tri) {
        glm::vec3 q = ClosestPointOnTriangle(p, tri.a, tri.b, tri.c);
        return glm::length(p - q);
    }

    bool MeshSDF::RayIntersectsTriangle(
        glm::vec3 const &    orig,
        glm::vec3 const &    dir,
        MeshTriangle const & tri) {
        float const eps = 1e-8f;

        glm::vec3 edge1 = tri.b - tri.a;
        glm::vec3 edge2 = tri.c - tri.a;
        glm::vec3 h     = glm::cross(dir, edge2);
        float     a     = glm::dot(edge1, h);

        if (std::abs(a) < eps)
            return false;

        float     f = 1.0f / a;
        glm::vec3 s = orig - tri.a;
        float     u = f * glm::dot(s, h);
        if (u < 0.0f || u > 1.0f)
            return false;

        glm::vec3 q = glm::cross(s, edge1);
        float     v = f * glm::dot(dir, q);
        if (v < 0.0f || u + v > 1.0f)
            return false;

        float t = f * glm::dot(edge2, q);
        return t > eps;
    }

    bool MeshSDF::IsInside(glm::vec3 const & p) const {
        glm::vec3 const dir(1.0f, 0.173f, 0.071f);
        int             intersections = 0;
        for (auto const & tri : _triangles) {
            if (RayIntersectsTriangle(p, dir, tri))
                ++intersections;
        }
        return (intersections & 1) == 1;
    }

    int MeshSDF::GridOffset(int x, int y, int z) const {
        return x + _gridResolution * (y + _gridResolution * z);
    }

    float MeshSDF::OutsideBoundsDistance(glm::vec3 const & localPoint) const {
        glm::vec3 const closest = glm::clamp(localPoint, _gridMin, _gridMax);
        return glm::length(localPoint - closest);
    }

    float MeshSDF::ExactSignedDistance(glm::vec3 const & localPoint) const {
        if (_triangles.empty())
            return std::numeric_limits<float>::max();

        float minDist = std::numeric_limits<float>::max();
        for (auto const & tri : _triangles)
            minDist = std::min(minDist, PointTriangleDistance(localPoint, tri));

        return IsInside(localPoint) ? -minDist : minDist;
    }

    void MeshSDF::BuildDistanceGrid() {
        _gridResolution = 32;
        _gridMin        = _minBounds - glm::vec3(0.08f);
        _gridMax        = _maxBounds + glm::vec3(0.08f);
        _gridDx         = (_gridMax - _gridMin) / float(_gridResolution - 1);
        _distanceGrid.assign(
            _gridResolution * _gridResolution * _gridResolution,
            std::numeric_limits<float>::max());

        for (int z = 0; z < _gridResolution; ++z) {
            for (int y = 0; y < _gridResolution; ++y) {
                for (int x = 0; x < _gridResolution; ++x) {
                    glm::vec3 const p =
                        _gridMin + glm::vec3(x, y, z) * _gridDx;
                    _distanceGrid[GridOffset(x, y, z)] = ExactSignedDistance(p);
                }
            }
        }
    }

    float MeshSDF::SignedDistance(glm::vec3 const & localPoint) const {
        if (_triangles.empty())
            return std::numeric_limits<float>::max();

        if (_distanceGrid.empty()
            || localPoint.x < _gridMin.x || localPoint.x > _gridMax.x
            || localPoint.y < _gridMin.y || localPoint.y > _gridMax.y
            || localPoint.z < _gridMin.z || localPoint.z > _gridMax.z) {
            return OutsideBoundsDistance(localPoint);
        }

        glm::vec3 const gridCoord =
            (localPoint - _gridMin) / glm::max(_gridDx, glm::vec3(1e-6f));
        glm::ivec3 const base = glm::clamp(
            glm::ivec3(glm::floor(gridCoord)),
            glm::ivec3(0),
            glm::ivec3(_gridResolution - 2));
        glm::vec3 const frac = glm::clamp(
            gridCoord - glm::vec3(base),
            glm::vec3(0.0f),
            glm::vec3(1.0f));

        auto sample = [&](int dx, int dy, int dz) {
            return _distanceGrid[GridOffset(base.x + dx, base.y + dy, base.z + dz)];
        };

        float const c00 = glm::mix(sample(0, 0, 0), sample(1, 0, 0), frac.x);
        float const c10 = glm::mix(sample(0, 1, 0), sample(1, 1, 0), frac.x);
        float const c01 = glm::mix(sample(0, 0, 1), sample(1, 0, 1), frac.x);
        float const c11 = glm::mix(sample(0, 1, 1), sample(1, 1, 1), frac.x);
        float const c0  = glm::mix(c00, c10, frac.y);
        float const c1  = glm::mix(c01, c11, frac.y);
        return glm::mix(c0, c1, frac.z);
    }

} // namespace VCX::Labs::Final
