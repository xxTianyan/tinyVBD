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
    std::vector<Vec3> p;          // 当前位置 (last frame pos)
    std::vector<Vec3> p_pred;     // 预测位置 (predict pos)
    std::vector<Vec3> p_inertia;  // 惯性预测 (inertia prediction)
    std::vector<Vec3> v;          // 速度
    std::vector<Vec3> accel;          // 加速度
    std::vector<Vec3> n;          // 法线
    std::vector<Vec3> inv_m;          // 质量

    [[nodiscard]] inline size_t size() const { return p.size(); }

    // 拓扑信息
    std::vector<tetrahedron> m_tets;
    std::vector<triangle> m_tris;
    std::vector<edge> m_edges;
    ForceElementAdjacencyInfo adjacencyInfo;

    // 渲染信息
    std::vector<VertexId> m_surface_tris;
    size_t base_offset = 0;

    // 重构后的 resize：从 18 行缩减到 6 行，极难出错
    void resize(const size_t n_nodes) {
        p.resize(n_nodes);
        p_pred.resize(n_nodes);
        p_inertia.resize(n_nodes);
        v.resize(n_nodes);
        accel.resize(n_nodes);
        n.resize(n_nodes);
        inv_m.resize(n_nodes);
    }

    // 清空数据
    void clear_topology() {
        m_edges.clear();
        m_tris.clear();
        m_tets.clear();
        m_surface_tris.clear();
    }
};

void BuildAdjacency(mesh_on_cpu& mesh);

void ParseMSH(const std::string& path, mesh_on_cpu* cpu_mesh);

std::vector<float> ComputeNormal(mesh_on_cpu* cpu_mesh);

IndexBuffer BuildSurfaceTriangles(const std::vector<tetrahedron>& tets);

IndexBuffer BuildSurfaceTriangles(const std::vector<triangle>& tris);

std::vector<float> assemble_vertices(const mesh_on_cpu* cpu_mesh);

#endif //MESH_H
