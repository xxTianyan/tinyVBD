//
// Created by 徐天焱 on 2025/11/4.
//

#ifndef MESHLOADER_HPP
#define MESHLOADER_HPP

#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include "Types.h"


struct mesh_on_cpu {
    std::vector<float> px, py, pz;  // last frame pos
    std::vector<float> px_pred, py_pred, pz_pred;  // predict pos
    std::vector<float> vx, vy, vz;
    std::vector<float> fx, fy, fz;
    std::vector<float> nx, ny, nz;  // normal

    [[nodiscard]] inline size_t size() const {return nx.size();}

    void resize(size_t n) {
        auto rsf = [&](auto& v){ v.resize(n); };
        rsf(px); rsf(py); rsf(pz);
        rsf(px_pred); rsf(py_pred); rsf(pz_pred);
        rsf(vx); rsf(vy); rsf(vz);
        rsf(fx); rsf(fy); rsf(fz);
        rsf(nx); rsf(ny); rsf(nz);
    }
};

struct mesh_on_gpu {
    std::vector<float> vertices;
    std::vector<float> norm;
};

struct tetrahedron {
    std::array<VertexId, 4> vertices{0,0,0,0};

    tetrahedron() = default;
    tetrahedron(const VertexId vtex1, const VertexId vtex2, const VertexId vtex3, const VertexId vtex4) :
    vertices{vtex1, vtex2, vtex3, vtex4} {};
};

struct NodeTetAdj {
    IndexBuffer offsets;
    IndexBuffer incidentTets;
};

NodeTetAdj buildNodeTetAdj(size_t num_nodes, const std::vector<tetrahedron>& tets);

std::vector<tetrahedron> ParseMSH(const std::string& path, mesh_on_cpu& cpu_mesh);

IndexBuffer ParseOBJ(const std::string& path, size_t num_nodes);

std::vector<float> ComputeNormal(mesh_on_cpu& cpu_mesh, const IndexBuffer& tri_indices);

IndexBuffer BuildSurfaceTriangles(const std::vector<tetrahedron>& tets);

std::vector<float> assemble_vertices(const mesh_on_cpu& cpu_mesh);

#endif //MESHLOADER_HPP
