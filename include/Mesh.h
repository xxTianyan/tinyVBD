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

struct mesh_on_cpu {
    std::vector<float> px, py, pz;  // last frame pos
    std::vector<float> px_pred, py_pred, pz_pred;  // predict pos
    std::vector<float> vx, vy, vz;
    std::vector<float> fx, fy, fz;
    std::vector<float> nx, ny, nz;  // normal

    // physical properties
    std::vector<float> inv_mass; // inverse of vertex mass
    std::vector<float> mass;     // vertex mass

    // material property / also can make global
    float mu = 1000.0f;    // Neo-Hookean paras
    float lambda = 4000.0f;
    bool inited = false;

    [[nodiscard]] inline size_t size() const {return nx.size();}

    std::vector<tetrahedron> m_tets_local;
    std::vector<VertexId> m_surface_tris_local;  // for rendering

    // pre-computing of tets
    std::vector<Mat3> tet_Dm_inv; // inverse of shape matrix Dm^-1
    std::vector<float> tet_vol;   // the volume of tet

    size_t base_offest = 0;

    void resize(size_t n) {
        auto rsf = [&](auto& v){ v.resize(n); };
        rsf(px); rsf(py); rsf(pz);
        rsf(px_pred); rsf(py_pred); rsf(pz_pred);
        rsf(vx); rsf(vy); rsf(vz);
        rsf(fx); rsf(fy); rsf(fz);
        rsf(nx); rsf(ny); rsf(nz);
        rsf(inv_mass); rsf(mass);
    }

    void InitializePhysics(float density);
};

struct NodeTetAdj {
    std::vector<uint32_t> offsets;        // n + 1 size, and tell which segment of incidentTets is relevant.
    std::vector<uint32_t> incidentTets;
};

inline bool isValidVertexId(const uint32_t value) {
    return value <= INVALID_VERTEX_ID && value > 0;  // > 0 , since index from msh file begins at 1
}

NodeTetAdj BuildNodeTetAdj(size_t num_nodes, const std::vector<tetrahedron>& tets);

void ParseMSH(const std::string& path, mesh_on_cpu* cpu_mesh);

std::vector<float> ComputeNormal(mesh_on_cpu* cpu_mesh);

IndexBuffer BuildSurfaceTriangles(const std::vector<tetrahedron>& tets);

std::vector<float> assemble_vertices(const mesh_on_cpu* cpu_mesh);

#endif //MESHLOADER_H
