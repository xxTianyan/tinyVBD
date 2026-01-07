//
// Created by tianyan on 12/22/25.
//

#include <stdexcept>
#include "VBDSolver.h"
#include "Math.hpp"

namespace {
    void AssignOffsets(const size_t num_nodes, std::vector<uint32_t>& offsets) {
        offsets.assign(num_nodes + 1, 0u);
    }

    template <class Elem, class GetVertex>
    void BuildVertexIncidentCSR(const size_t num_nodes,
                                const std::vector<Elem>& elems,
                                const uint32_t verts_per_elem,
                                GetVertex get_vertex,
                                AdjacencyCSR& adj) {
        auto& offsets = adj.offsets;
        auto& incidents = adj.incidents;

        offsets.assign(num_nodes + 1, 0u);
        incidents.clear();

        if (elems.empty()) {
            return;
        }

        if (verts_per_elem > 4u) {
            throw std::runtime_error("verts_per_elem > 4 is not supported by AdjacencyCSR::pack");
        }

        for (uint32_t elem_id = 0; elem_id < static_cast<uint32_t>(elems.size()); ++elem_id) {
            const auto& elem = elems[elem_id];
            for (uint32_t k = 0; k < verts_per_elem; ++k) {
                const auto v = static_cast<uint32_t>(get_vertex(elem, k));
                offsets[v + 1] += 1u;
            }
        }

        for (size_t i = 1; i < offsets.size(); ++i) {
            offsets[i] += offsets[i - 1];
        }

        incidents.resize(offsets.back());
        std::vector<uint32_t> cursor = offsets;

        for (uint32_t elem_id = 0; elem_id < static_cast<uint32_t>(elems.size()); ++elem_id) {
            const auto& elem = elems[elem_id];
            for (uint32_t k = 0; k < verts_per_elem; ++k) {
                const auto v = static_cast<uint32_t>(get_vertex(elem, k));
                const uint32_t dst = cursor[v]++;
                incidents[dst] = AdjacencyCSR::pack(elem_id, k);
            }
        }
    }

    void BuildAdjacency(const MModel& model, ForceElementAdjacencyInfo& adjacency_info) {
        const size_t num_nodes = model.total_particles();
        if (!model.edges.empty()) {
            BuildVertexIncidentCSR(
                num_nodes,
                model.edges,
                4u,
                [](const edge& e, uint32_t k) { return static_cast<uint32_t>(e.vertices[k]); },
                adjacency_info.vertex_edges
            );
        } else {
            AssignOffsets(num_nodes, adjacency_info.vertex_edges.offsets);
            adjacency_info.vertex_edges.incidents.clear();
        }

        if (!model.tris.empty()) {
            BuildVertexIncidentCSR(
                num_nodes,
                model.tris,
                3u,
                [](const triangle& t, uint32_t k) { return static_cast<uint32_t>(t.vertices[k]); },
                adjacency_info.vertex_faces
            );
        } else {
            AssignOffsets(num_nodes, adjacency_info.vertex_faces.offsets);
            adjacency_info.vertex_faces.incidents.clear();
        }

        if (!model.tets.empty()) {
            BuildVertexIncidentCSR(
                num_nodes,
                model.tets,
                4u,
                [](const tetrahedron& t, uint32_t k) { return static_cast<uint32_t>(t.vertices[k]); },
                adjacency_info.vertex_tets
            );
        } else {
            AssignOffsets(num_nodes, adjacency_info.vertex_tets.offsets);
            adjacency_info.vertex_tets.incidents.clear();
        }
    }
}

void VBDSolver::Init(const Scene& scene) {
    const auto& model = scene.model_;
    const size_t num_nodes = model.total_particles();
    inertia_.assign(num_nodes, Vec3::Zero());
    prev_pos_.assign(num_nodes, Vec3::Zero());

    if (model.topology_version != topology_version_
        || adjacency_info_.vertex_faces.offsets.size() != num_nodes + 1) {
        BuildAdjacency(model, adjacency_info_);
        topology_version_ = model.topology_version;
    }
}

void VBDSolver::Step(Scene& scene, const float dt) {
    if (dt <= 0.0f) {
        return;
    }

    Init(scene);

    forward_step(scene, dt);
    for (int iter = 0; iter < num_iters; ++iter) {
        solve(scene, dt);
    }
    update_velocity(scene, dt);
    scene.state_in_ = scene.state_out_;
}

void VBDSolver::accumulate_stvk_triangle_force_hessian(const std::span<const Vec3> pos,
                                                   const MMaterial& mat,
                                                   const triangle& face,
                                                   const uint32_t vtex_order,
                                                   Vec3& force,
                                                   Mat3& H) {
    // advised by newton physics, evaluate_stvk_force_hessian function
    // StVK energy density: psi = mu * ||G||_F^2 + 0.5 * lambda * (trace(G))^2

    if (vtex_order > 2)
        throw std::runtime_error("vtex order is over stvk triangle limt");

    const Vec3& x0 = pos[face.vertices[0]];
    const Vec3& x1 = pos[face.vertices[1]];
    const Vec3& x2 = pos[face.vertices[2]];
    const auto mu = mat.mu();
    const auto lambda = mat.lambda();

    // Deformation gradient F = [f0, f1] (3x2 matrix as two 3D column vectors)
    const auto DmInv00 = face.Dm_inv(0,0);
    const auto DmInv01 = face.Dm_inv(0,1);
    const auto DmInv10 = face.Dm_inv(1,0);
    const auto DmInv11 = face.Dm_inv(1,1);

    // Compute F columns directly: F = [x01, x02] * tri_pose = [f0, f1]
    const Vec3 f0 = (x1 - x0) * DmInv00 + (x2 - x0) * DmInv10;
    const Vec3 f1 = (x1 - x0) * DmInv01 + (x2 - x0) * DmInv11;

    // Green strain tensor: G = 0.5(F^T F - I) = [[G00, G01], [G01, G11]] (symmetric 2x2)
    const auto f0_dot_f0 = f0.dot(f0);
    const auto f1_dot_f1 = f1.dot(f1);
    const auto f0_dot_f1 = f0.dot(f1);

    const auto G00 = 0.5 * (f0_dot_f0 - 1.0);
    const auto G11 = 0.5 * (f1_dot_f1 - 1.0);
    const auto G01 = 0.5 * f0_dot_f1;

    // Frobenius norm squared of Green strain: ||G||_F^2 = G00^2 + G11^2 + 2 * G01^2
    float G_frobenius_sq = G00 * G00 + G11 * G11 + 2.0f * G01 * G01;
    if (G_frobenius_sq < 1.0e-20) {
        return;
    }

    const float trace_G = G00 + G11;

    // First Piola-Kirchhoff stress tensor (StVK model)
    // PK1 = 2*mu*F*G + lambda*trace(G)*F = [PK1_col0, PK1_col1] (3x2)
    const auto lambda_trace_G = lambda * trace_G;
    const auto two_mu = 2.0f * mu;

    const Vec3 PK1_col0 = f0 * (two_mu * G00 + lambda_trace_G) + f1 * (two_mu * G01);
    const Vec3 PK1_col1 = f0 * (two_mu * G01) + f1 * (two_mu * G11 + lambda_trace_G);

    const auto mask0 = static_cast<float>(vtex_order == 0);
    const auto mask1 = static_cast<float>(vtex_order == 1);
    const auto mask2 = static_cast<float>(vtex_order == 2);

    // Deformation gradient derivatives w.r.t. current vertex position
    const auto df0_dx = DmInv00 * (mask1 - mask0) + DmInv10 * (mask2 - mask0);
    const auto df1_dx = DmInv01 * (mask1 - mask0) + DmInv11 * (mask2 - mask0);

    // Force via chain rule: force = -(dpsi/dF) : (dF/dx)
    const Vec3 dpsi_dx = PK1_col0 * df0_dx + PK1_col1 * df1_dx;
    Vec3 delta_force = -dpsi_dx;

    // Hessian computation using Cauchy-Green invariants
    const auto df0_dx_sq = df0_dx * df0_dx;
    const auto df1_dx_sq = df1_dx * df1_dx;
    const auto  df0_df1_cross = df0_dx * df1_dx;

    const auto Ic = f0_dot_f0 + f1_dot_f1;
    const auto two_dpsi_dIc = -mu + (0.5 * Ic - 1.0) * lambda;
    const Mat3 I33 = Mat3::Identity();

   const Mat3 f0_outer_f0 = f0 * f0.transpose();
   const Mat3 f1_outer_f1 = f1 * f1.transpose();
   const Mat3 f0_outer_f1 = f0 * f1.transpose();
   const Mat3 f1_outer_f0 = f1 * f0.transpose();

   const Mat3 H_IIc00_scaled = mu * (f0_dot_f0 * I33 + 2.0 * f0_outer_f0 + f1_outer_f1);
   const Mat3 H_IIc11_scaled = mu * (f1_dot_f1 * I33 + 2.0 * f1_outer_f1 + f0_outer_f0);
   const Mat3 H_IIc01_scaled = mu * (f0_dot_f1 * I33 + f1_outer_f0);

    // d2(psi)/dF^2 components
   const Mat3 d2E_dF2_00 = lambda * f0_outer_f0 + two_dpsi_dIc * I33 + H_IIc00_scaled;
   const Mat3 d2E_dF2_01 = lambda * f0_outer_f1 + H_IIc01_scaled;
   const Mat3 d2E_dF2_11 = lambda * f1_outer_f1 + two_dpsi_dIc * I33 + H_IIc11_scaled;

    // Chain rule: H = (dF/dx)^T * (d2(psi)/dF^2) * (dF/dx)
   Mat3 delta_hessian = df0_dx_sq * d2E_dF2_00 + df1_dx_sq * d2E_dF2_11 +
       df0_df1_cross * (d2E_dF2_01 + d2E_dF2_01.transpose());

    force += delta_force * face.rest_area;
    H += delta_hessian * face.rest_area;
}

void VBDSolver::accumulate_dihedral_angle_based_bending_force_hessian(const std::span<const Vec3> pos,
                                                                const MMaterial& mat,
                                                                const edge& e,
                                                                const uint32_t vtex_order,
                                                                Vec3& force,
                                                                Mat3& H) {
    // advised by function with the same name in newton physics.
    constexpr float eps = 1e-6f;

    const auto& x0 = pos[e.vertices[0]];
    const auto& x1 = pos[e.vertices[1]];
    const auto& x2 = pos[e.vertices[2]];
    const auto& x3 = pos[e.vertices[3]];

    // Compute edge vectors
    const Vec3 x02 = x2 - x0;
    const Vec3 x03 = x3 - x0;
    const Vec3 x12 = x2 - x1;
    const Vec3 x13 = x3 - x1;
    const Vec3 edge_line = x3 - x2;

    // Compute normals
    const Vec3 n1 = x02.cross(x03);
    const Vec3 n2 = x13.cross(x12);

    const float n1_norm = n1.norm();
    const float n2_norm = n2.norm();
    const float e_norm = edge_line.norm();

    if (n1_norm < eps || n2_norm < eps) return;

    const Vec3 n1_hat = n1 / n1_norm;
    const Vec3 n2_hat = n2 / n2_norm;
    const Vec3 e_hat = edge_line / e_norm;

    const auto sin_theta = (n1_hat.cross(n2_hat)).dot(e_hat);
    const auto cos_theta = n1_hat.dot(n2_hat);

    const auto theta = std::atan2(sin_theta, cos_theta);
    const auto k = mat.bend_stiff() * e.rest_length;
    const auto dE_dtheta = k * (theta - e.rest_theta);

    // Pre-compute skew matrices (shared across all angle derivative computations)
    Mat3 skew_e = TY::Skew(edge_line);
    Mat3 skew_x02 = TY::Skew(x02);
    Mat3 skew_x03 = TY::Skew(x03);
    Mat3 skew_x12 = TY::Skew(x12);
    Mat3 skew_x13 = TY::Skew(x13);
    Mat3 skew_n1 = TY::Skew(n1_hat);
    Mat3 skew_n2 = TY::Skew(n2_hat);

    // Compute the derivatives of unit normals with respect to each vertex; required for computing angle derivatives
    const auto dn1hat_dx0 = TY::ComputeNormalizedVectorDerivative(n1_norm, n1_hat, skew_e);
    const Mat3 dn2hat_dx0 = Mat3::Zero();

    const Mat3 dn1hat_dx1 = Mat3::Zero();
    const Mat3 dn2hat_dx1 = TY::ComputeNormalizedVectorDerivative(n2_norm, n2_hat, -skew_e);

    const Mat3 dn1hat_dx2 = TY::ComputeNormalizedVectorDerivative(n1_norm, n1_hat, -skew_x03);
    const Mat3 dn2hat_dx2 = TY::ComputeNormalizedVectorDerivative(n2_norm, n2_hat, skew_x13);

    const Mat3 dn1hat_dx3 = TY::ComputeNormalizedVectorDerivative(n1_norm, n1_hat, skew_x02);
    const Mat3 dn2hat_dx3 = TY::ComputeNormalizedVectorDerivative(n2_norm, n2_hat, -skew_x12);

    // Compute all angle derivatives (required for damping)
    const Vec3 dtheta_dx0 = TY::ComputeAngleDerivative(
    n1_hat, n2_hat, e_hat, dn1hat_dx0, dn2hat_dx0, sin_theta, cos_theta, skew_n1, skew_n2
    );

    const Vec3 dtheta_dx1 = TY::ComputeAngleDerivative(
    n1_hat, n2_hat, e_hat, dn1hat_dx1, dn2hat_dx1, sin_theta, cos_theta, skew_n1, skew_n2
    );

    const Vec3 dtheta_dx2 = TY::ComputeAngleDerivative(
    n1_hat, n2_hat, e_hat, dn1hat_dx2, dn2hat_dx2, sin_theta, cos_theta, skew_n1, skew_n2
    );

    const Vec3 dtheta_dx3 = TY::ComputeAngleDerivative(
    n1_hat, n2_hat, e_hat, dn1hat_dx3, dn2hat_dx3, sin_theta, cos_theta, skew_n1, skew_n2
    );

    // Use float masks for branch-free selection
    const auto mask0 = static_cast<float>(vtex_order == 0);
    const auto mask1 = static_cast<float>(vtex_order == 1);
    const auto mask2 = static_cast<float>(vtex_order == 2);
    const auto mask3 = static_cast<float>(vtex_order == 3);

    // Select the derivative for the current vertex without branching
    const Vec3 dtheta_dx = dtheta_dx0 * mask0 + dtheta_dx1 * mask1 + dtheta_dx2 * mask2 + dtheta_dx3 * mask3;

    // Compute elastic force and hessian
    const Vec3 bending_force = -dE_dtheta * dtheta_dx;
    const Mat3 bending_hessian = k * dtheta_dx * dtheta_dx.transpose();

    force += bending_force;
    H += bending_hessian;
}

void VBDSolver::forward_step(Scene& scene, const float dt) {
    const auto& model = scene.model_;
    const auto& state_in = scene.state_in_;
    auto& state_out = scene.state_out_;

    const size_t num_nodes = model.total_particles();
    if (state_out.particle_pos.size() != num_nodes) {
        state_out.resize(num_nodes);
    }

    state_out.particle_pos = state_in.particle_pos;
    state_out.particle_vel = state_in.particle_vel;
    state_out.particle_force = state_in.particle_force;

    for (size_t i = 0; i < num_nodes; ++i) {
        prev_pos_[i] = state_in.particle_pos[i];
        const float inv_mass = model.particle_inv_mass[i];
        Vec3 accel = Vec3::Zero();
        if (inv_mass > 0.0f) {
            accel = state_in.particle_force[i] * inv_mass;
        }
        inertia_[i] = state_in.particle_pos[i] + dt * state_in.particle_vel[i] + dt * dt * accel;
    }
}

void VBDSolver::solve(Scene& scene, const float dt) {
    const auto& model = scene.model_;
    auto& state_out = scene.state_out_;
    const VertexID num_nodes = static_cast<VertexID>(model.total_particles());

    for (VertexID vtex_id = 0; vtex_id < num_nodes; ++vtex_id) {
        auto& pos = state_out.particle_pos[vtex_id];
        const auto& inv_mass = model.particle_inv_mass[vtex_id];

        if (inv_mass <= 0.0f) {
            continue;
        }

        const auto& face_adjacency = adjacency_info_.vertex_faces;
        const auto& edge_adjacency = adjacency_info_.vertex_edges;

        Vec3 force = Vec3::Zero();
        Mat3 hessian = Mat3::Zero();

        force += -(pos - inertia_[vtex_id]) / (inv_mass * dt * dt);
        hessian += Mat3::Identity() / (inv_mass * dt * dt);

        for (uint32_t f = face_adjacency.begin(vtex_id); f < face_adjacency.end(vtex_id); ++f) {
            const auto pack = face_adjacency.incidents[f];
            const auto face_id = AdjacencyCSR::unpack_id(pack);
            const auto order = AdjacencyCSR::unpack_order(pack);
            const auto& face = model.tris[face_id];
            accumulate_stvk_triangle_force_hessian(state_out.particle_pos, material_, face, order, force, hessian);
        }

        for (uint32_t e = edge_adjacency.begin(vtex_id); e < edge_adjacency.end(vtex_id); ++e) {
            const auto pack = edge_adjacency.incidents[e];
            const auto edge_id = AdjacencyCSR::unpack_id(pack);
            const auto order = AdjacencyCSR::unpack_order(pack);
            const auto& edge = model.edges[edge_id];
            accumulate_dihedral_angle_based_bending_force_hessian(state_out.particle_pos, material_, edge, order, force, hessian);
        }

        const auto delta_x = TY::SolveSPDOrRegularize(hessian, force);
        pos += delta_x;
    }
}

void VBDSolver::update_velocity(Scene& scene, const float dt) {
    const auto& model = scene.model_;
    auto& state_out = scene.state_out_;
    const size_t num_nodes = model.total_particles();

    for (size_t i = 0; i < num_nodes; ++i) {
        const float inv_mass = model.particle_inv_mass[i];
        if (inv_mass <= 0.0f) {
            state_out.particle_vel[i] = Vec3::Zero();
            continue;
        }
        state_out.particle_vel[i] = (state_out.particle_pos[i] - prev_pos_[i]) / dt;
    }
}

























