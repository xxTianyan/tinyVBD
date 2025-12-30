//
// Created by 徐天焱 on 2025/11/4.
//

#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include <array>
#include "Types.h"
#include "AdjacencyCSR.hpp"

struct tetrahedron {
    std::array<VertexId, 4> vertices{0,0,0,0};
    tetrahedron(const VertexId vtex0, const VertexId vtex1, const VertexId vtex2, const VertexId vtex3,
        const Vec3& vtex0_pos, const Vec3& vtex1_pos, const Vec3& vtex2_pos, const Vec3& vtex3_pos) :
    vertices{vtex0, vtex1, vtex2, vtex3} {


    };
};

struct triangle {
    std::array<VertexId, 3> vertices{0,0,0};
    float rest_area;
    Mat2 Dm_inv;

    triangle(const VertexId vtex0, const VertexId vtex1, const VertexId vtex2, const Vec3& vtex0_pos, const Vec3& vtex1_pos, const Vec3& vtex2_pos) :
    vertices{vtex0, vtex1, vtex2} {

    };
};

struct edge {
    std::array<VertexId, 4> vertices{0,0,0,0};
    edge(const VertexId vtex0, const VertexId vtex1, const VertexId vtex2, const VertexId vtex3,
        const Vec3& vtex0_pos, const Vec3& vtex1_pos, const Vec3& vtex2_pos, const Vec3& vtex3_pos) :
    vertices{vtex0, vtex1, vtex2, vtex3} {

    };

};

struct mesh_on_cpu {
    std::vector<Vec3> pos;          // 当前位置 (current frame pos)
    std::vector<Vec3> prev_pos;     // 预测位置 (last frame pos)
    std::vector<Vec3> inertia_pos;  // 惯性预测 (inertia prediction)
    std::vector<Vec3> vel;          // 速度
    std::vector<Vec3> accel;          // 加速度
    std::vector<Vec3> n;          // 法线
    std::vector<float> inv_mass;          // 质量

    [[nodiscard]] inline size_t size() const { return pos.size(); }

    // 拓扑信息
    std::vector<tetrahedron> m_tets;
    std::vector<triangle> m_tris;
    std::vector<edge> m_edges;
    ForceElementAdjacencyInfo adjacencyInfo;

    // 渲染信息
    std::vector<VertexId> m_surface_tris;
    size_t base_offset = 0;

    // resize
    void resize(const size_t n_nodes) {
        pos.resize(n_nodes);
        prev_pos.resize(n_nodes);
        inertia_pos.resize(n_nodes);
        vel.resize(n_nodes);
        accel.resize(n_nodes);
        n.resize(n_nodes);
        inv_mass.resize(n_nodes);
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

void DistributeMass(mesh_on_cpu& mesh);

void InitMesh(mesh_on_cpu& mesh);

void ParseMSH(const std::string& path, mesh_on_cpu* cpu_mesh);

std::vector<float> ComputeNormal(mesh_on_cpu* cpu_mesh);

IndexBuffer BuildSurfaceTriangles(const std::vector<tetrahedron>& tets);

IndexBuffer BuildSurfaceTriangles(const std::vector<triangle>& tris);

std::vector<float> AssembleVertices(const mesh_on_cpu* cpu_mesh);

#endif //MESH_H
