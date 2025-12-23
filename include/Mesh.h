//
// Created by 徐天焱 on 2025/11/4.
//

#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include <array>
#include "Types.h"

struct tetrahedron {
    std::array<VertexId, 4> vertices{0,0,0,0};
    tetrahedron() = default;
    tetrahedron(const VertexId vtex1, const VertexId vtex2, const VertexId vtex3, const VertexId vtex4) :
    vertices{vtex1, vtex2, vtex3, vtex4} {};
};

struct triangle {
    std::array<VertexId, 3> vertices{0,0,0};
    triangle() = default;
    triangle(const VertexId vtex1, const VertexId vtex2, const VertexId vtex3) :
    vertices{vtex1, vtex2, vtex3} {};
};

struct edge {
    std::array<VertexId, 2> vertices{0,0};
    edge() = default;
    edge(const VertexId vtex1, const VertexId vtex2) :
    vertices{vtex1, vtex2} {};

};

struct ForceElementAdjacencyInfo{
    std::vector<uint32_t> v_adj_edges_offsets;
    std::vector<uint32_t> v_adj_edges;

    std::vector<uint32_t> v_adj_faces_offsets;
    std::vector<uint32_t> v_adj_faces;

    std::vector<uint32_t> v_adj_tet_offsets;
    std::vector<uint32_t> v_adj_tets;
};

struct mesh_on_cpu {
    std::vector<float> px, py, pz;  // last frame pos
    std::vector<float> px_pred, py_pred, pz_pred;  // predict pos
    std::vector<float> inertia_x, inertia_y, inertia_z;  // inertia prediction
    std::vector<float> vx, vy, vz;
    std::vector<float> fx, fy, fz;
    std::vector<float> nx, ny, nz;  // normal
    Dimension dim;

    [[nodiscard]] inline size_t size() const {return nx.size();}

    // for physical simulations
    std::vector<tetrahedron> m_tets;  // empty if it's 2d mesh
    std::vector<triangle> m_tris;
    std::vector<edge> m_edges;
    ForceElementAdjacencyInfo adjacencyInfo;

    // for rendering
    std::vector<uint32_t> m_surface_tris;

    size_t base_offest = 0;

    void resize(size_t n) {
        auto rsf = [&](auto& v){ v.resize(n); };
        rsf(px); rsf(py); rsf(pz);
        rsf(px_pred); rsf(py_pred); rsf(pz_pred);
        rsf(inertia_x); rsf(inertia_y); rsf(inertia_z);
        rsf(vx); rsf(vy); rsf(vz);
        rsf(fx); rsf(fy); rsf(fz);
        rsf(nx); rsf(ny); rsf(nz);
    }

};

inline bool isValidVertexId(const uint32_t value) {
    return value <= INVALID_VERTEX_ID && value > 0;  // > 0 , since index from msh file begins at 1
}

void BuildAdjacency(mesh_on_cpu* mesh);

void ParseMSH(const std::string& path, mesh_on_cpu* cpu_mesh);

std::vector<float> ComputeNormal(mesh_on_cpu* cpu_mesh);

IndexBuffer BuildSurfaceTriangles(const std::vector<tetrahedron>& tets);

IndexBuffer BuildSurfaceTriangles(const std::vector<triangle>& tris);

std::vector<float> assemble_vertices(const mesh_on_cpu* cpu_mesh);

#endif //MESH_H
