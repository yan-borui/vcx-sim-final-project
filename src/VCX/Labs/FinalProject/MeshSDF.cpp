#include "Labs/FinalProject/MeshSDF.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#include <tiny_obj_loader.h>

namespace VCX::Labs::Final {

    namespace {
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
                }

                _triangles.push_back({ verts[0], verts[1], verts[2] });
                indexOffset += 3;
            }
        }

        std::cout << "MeshSDF loaded triangles: " << _triangles.size()
                  << " from " << filePath << std::endl;

        return ! _triangles.empty();
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

    float MeshSDF::SignedDistance(glm::vec3 const & localPoint) const {
        if (_triangles.empty())
            return std::numeric_limits<float>::max();

        float minDist = std::numeric_limits<float>::max();
        for (auto const & tri : _triangles)
            minDist = std::min(minDist, PointTriangleDistance(localPoint, tri));

        return IsInside(localPoint) ? -minDist : minDist;
    }

} // namespace VCX::Labs::Final