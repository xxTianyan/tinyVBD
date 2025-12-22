//
// Created by 徐天焱 on 2025/11/4.
//

#ifndef MESHLOADER_H
#define MESHLOADER_H

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

    // physical properties
    std::vector<float> inv_mass; // inverse of vertex mass
    std::vector<float> mass;     // vertex mass

    [[nodiscard]] inline size_t size() const {return nx.size();}

    // for physical simulations
    std::vector<tetrahedron> m_tets_local;  // empty if it's 2d mesh
    std::vector<triangle> m_tets;
    std::vector<edge> m_edges;

    // for rendering
    std::vector<VertexId> m_surface_tris_local;

    // pre-computing of tets
    std::vector<Mat3> tet_Dm_inv; // inverse of shape matrix Dm^-1
    std::vector<float> tet_vol;   // the volume of tet

    size_t base_offest = 0;

    void resize(size_t n) {
        auto rsf = [&](auto& v){ v.resize(n); };
        rsf(px); rsf(py); rsf(pz);
        rsf(px_pred); rsf(py_pred); rsf(pz_pred);
        rsf(inertia_x); rsf(inertia_y); rsf(inertia_z);
        rsf(vx); rsf(vy); rsf(vz);
        rsf(fx); rsf(fy); rsf(fz);
        rsf(nx); rsf(ny); rsf(nz);
        rsf(inv_mass); rsf(mass);
    }

    void InitializePhysics(float density);
};

inline bool isValidVertexId(const uint32_t value) {
    return value <= INVALID_VERTEX_ID && value > 0;  // > 0 , since index from msh file begins at 1
}

void BuildAdjacency(size_t num_nodes, const std::vector<tetrahedron>& tets);

void ParseMSH(const std::string& path, mesh_on_cpu* cpu_mesh);

std::vector<float> ComputeNormal(mesh_on_cpu* cpu_mesh);

IndexBuffer BuildSurfaceTriangles(const std::vector<tetrahedron>& tets);

std::vector<float> assemble_vertices(const mesh_on_cpu* cpu_mesh);

#endif //MESHLOADER_H
