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
        ;
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

// currently, state and model only store one deformable  mesh

struct State {
    std::vector<Vec3> particle_pos;          // current frame pos
    std::vector<Vec3> particle_vel;          // velocity
    std::vector<Vec3> particle_force;          // force

    void resize(const size_t n_nodes) {
        particle_pos.resize(n_nodes);
        particle_vel.resize(n_nodes);
        particle_force.resize(n_nodes);
    }

    void clear() {
        particle_pos.clear();
        particle_vel.clear();
        particle_force.clear();
    }

    void clean_force() {
        std::fill(particle_force.begin(), particle_force.end(), Vec3{0.f,0.0f,0.0f});
    }
};

struct MModel {
    std::vector<tetrahedron> tets;
    std::vector<triangle> tris;
    std::vector<edge> edges;
    std::vector<float> particle_inv_mass;    // inverse mass
    std::vector<uint8_t> if_particle_fixed;     // if fixed
    std::vector<VertexID> surface_tris;     // surface index buffer (triangles), size must be multiple of 3
    std::vector<Vec3> pos0; // initial positions
    std::vector<Vec3> vel0; // initial velocities (optional; default zero)

    size_t num_nodes = 0;      // decided when model is finalized
    size_t base_offset = 0;    // assigned by Scene when packing into global arrays (optional)
    [[nodiscard]] inline size_t size() const { return num_nodes; }

    // topology has no default constructor and size is uncertain, thus not reserve space here
    void init(const size_t n_nodes) {
        num_nodes = n_nodes;
        particle_inv_mass.resize(n_nodes);
        if_particle_fixed.resize(n_nodes);
        pos0.resize(n_nodes);
        vel0.resize(n_nodes);

        // 给出稳妥默认值（builder 可覆盖）
        std::fill(particle_inv_mass.begin(), particle_inv_mass.end(), 1.0f);
        std::fill(if_particle_fixed.begin(), if_particle_fixed.end(), uint8_t{0});
        std::fill(pos0.begin(), pos0.end(), Vec3::Zero());
        std::fill(vel0.begin(), vel0.end(), Vec3::Zero());
    }

    void clear_topology() {
        edges.clear();
        tris.clear();
        tets.clear();
        surface_tris.clear();
    }
};



// void ParseMSH(const std::string& path, mesh_on_cpu* cpu_mesh);

// IndexBuffer BuildSurfaceTriangles(const std::vector<tetrahedron>& tets);

// IndexBuffer BuildSurfaceTriangles(const std::vector<triangle>& tris);



#endif //TAIYI_MESH_H
