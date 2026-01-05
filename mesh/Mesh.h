//
// Created by 徐天焱 on 2025/11/4.
//

#ifndef TAIYI_MESH_H
#define TAIYI_MESH_H

#include <vector>
#include <string>
#include <array>
#include "Types.h"
#include "AdjacencyCSR.hpp"

struct tetrahedron {
    std::array<VertexID, 4> vertices{0,0,0,0};
    tetrahedron(const VertexID vtex0, const VertexID vtex1, const VertexID vtex2, const VertexID vtex3,
        const Vec3& vtex0_pos, const Vec3& vtex1_pos, const Vec3& vtex2_pos, const Vec3& vtex3_pos) :
    vertices{vtex0, vtex1, vtex2, vtex3} {


    };
};

struct triangle {
    std::array<VertexID, 3> vertices{0,0,0};
    float rest_area;
    Mat2 Dm_inv;

    triangle(const VertexID vtex0, const VertexID vtex1, const VertexID vtex2, const Vec3& vtex0_pos, const Vec3& vtex1_pos, const Vec3& vtex2_pos) :
    vertices{vtex0, vtex1, vtex2} {

        /*
         * advised by newton add_triangle function in builder.py
         */

        // build local 2d material basis from rest position
        const Vec3 e01 = vtex1_pos - vtex0_pos;
        const float L = e01.norm();
        constexpr float eps = 1e-8f;

        if (L < eps) {
            // Degenerate: vtex0_pos == vtex1_pos
            rest_area = 0.0f;
            Dm_inv.setZero();
            return;
        }

        const Vec3 e1 = e01 / L;

        const Vec3 e02 = vtex2_pos - vtex0_pos;
        Vec3 n = e01.cross(e02);
        const float n_len = n.norm();
        if (n_len < eps) {
            // Degenerate: collinear triangle in rest pose
            rest_area = 0.0f;
            Dm_inv.setZero();
            return;
        }

        n /= n_len;
        const Vec3 e2 = n.cross(e1);  // orthonormal within triangle plane

        // --- Rest 2D coordinates ---
        // u0 = (0,0)
        // u1 = (L,0)
        // u2 = (u,v)
        const float u = e02.dot(e1);
        const float v = e02.dot(e2);

        // Dm = [[L, u],
        //       [0, v]]
        // det(Dm) = L*v
        const float det = L * v;
        if (det < eps) {
            rest_area = 0.0f;
            Dm_inv.setZero();
            return;
        }

        // --- rest area ---
        rest_area = 0.5f * det;


        // --- explicit inverse of upper-triangular Dm ---
        // Dm^{-1} = [[ 1/L,   -u/(L*v)],
        //           [ 0,      1/v      ]]
        Dm_inv(0, 0) = 1.0f / L;
        Dm_inv(0, 1) = -u / det;  // -u/(L*v)
        Dm_inv(1, 0) = 0.0f;
        Dm_inv(1, 1) = 1.0f / v;
    };
};

struct edge {
    // convention [opp0, opp1, edge_start, edge_end]
    std::array<VertexID, 4> vertices{0,0,0,0};
    float rest_theta;
    float rest_length;
    edge(const VertexID vtex_opp0, const VertexID vtex_opp1, const VertexID vtex_e0, const VertexID vtex_e1,
        const Vec3& vtex_opp0_pos, const Vec3& vtex_opp1_pos, const Vec3& vtex_e0_pos, const Vec3& vtex_e1_pos) :
    vertices{vtex_opp0, vtex_opp1, vtex_e0, vtex_e1} {
        rest_length = (vtex_e1_pos - vtex_e0_pos).norm();
        rest_theta = ComputeRestDihedralAngle(vtex_opp0_pos, vtex_opp1_pos, vtex_e0_pos, vtex_e1_pos);
    };

    // the same definition with newton warp atan2(sin, cos)
    static float ComputeRestDihedralAngle(const Vec3& x0, const Vec3& x1, const Vec3& x2, const Vec3& x3) {

        constexpr float eps = 1e-6f;
        const Vec3 e = x3 - x2;
        const float e_norm = e.norm();
        if (e_norm < eps) return 0.0f;

        const Vec3 x02 = x2 - x0;
        const Vec3 x03 = x3 - x0;
        const Vec3 x12 = x2 - x1;
        const Vec3 x13 = x3 - x1;

        const Vec3 n1 = x02.cross(x03);
        const Vec3 n2 = x13.cross(x12);

        const float n1_norm = n1.norm();
        const float n2_norm = n2.norm();

        if (n1_norm < eps || n2_norm < eps) return 0.0f;

        const Vec3 n1_hat = n1 / n1_norm;
        const Vec3 n2_hat = n2 / n2_norm;
        const Vec3 e_hat = e / e_norm;

        const float sin_theta = (n1_hat.cross(n2_hat)).dot(e_hat);
        const float cos_theta = n1_hat.dot(n2_hat);

        return std::atan2(sin_theta, cos_theta);
    };
};

struct mesh_on_cpu {
    std::vector<Vec3> pos;          // current frame pos
    std::vector<Vec3> prev_pos;     // last frame pos
    std::vector<Vec3> inertia_pos;  // inertia prediction
    std::vector<Vec3> vel;          // velocity
    std::vector<Vec3> accel;          // acceleration
    std::vector<Vec3> n;          // normal
    std::vector<float> inv_mass;          // inverse mass

    std::vector<uint8_t> fixed;  // if fixed

    [[nodiscard]] inline size_t size() const { return pos.size(); }

    // topology information
    std::vector<tetrahedron> m_tets;
    std::vector<triangle> m_tris;
    std::vector<edge> m_edges;
    ForceElementAdjacencyInfo adjacencyInfo;  // maybe it's better to make it a unique_ptr

    // rendering information
    std::vector<VertexID> m_surface_tris;
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
        fixed.resize(n_nodes);
    }

    // clear topology
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

#endif //TAIYI_MESH_H
