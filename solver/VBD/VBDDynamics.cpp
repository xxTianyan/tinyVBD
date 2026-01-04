//
// Created by tianyan on 12/22/25.
//

#include "VBDDynamics.h"
#include <iostream>

#include "Math.hpp"

static void accumulate_stvk_triangle_force_hessian(const std::span<const Vec3> pos,
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

static void accumulate_dihedral_angle_based_bending_force_hessian(const std::span<const Vec3> pos,
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


void VBDSolver::forward_step(SimView &view, const float dt) {
    const size_t num_nodes = view.pos.size();
    for (size_t i = 0; i < num_nodes; ++i) {
        view.prev_pos[i] = view.pos[i];  // record previous pos
        view.inertia_pos[i] = view.pos[i] + dt * view.vel[i] + dt * dt * (view.accel[i] + view.gravity);
    }
}

void VBDSolver::solve(SimView& view, const float dt) {
    const VertexID num_nodes = view.pos.size();  // the number of nodes always should less than 0xFFFF

    // solve
    for (VertexID vtex_id = 0; vtex_id < num_nodes; ++vtex_id) {

        if (view.fixed[vtex_id]) continue;  // fixed flag up

        auto& pos = view.pos[vtex_id];
        const auto& inv_mass = view.inv_mass[vtex_id];

        if (inv_mass <= 0.0f) continue;

        const auto& inertia_pos = view.inertia_pos[vtex_id];
        const auto& prev_pos = view.prev_pos[vtex_id];
        // const auto& edge_adjacency = view.adj.vertex_edges;
        const auto& face_adjacency = view.adj.vertex_faces;
        const auto& edge_adjacency = view.adj.vertex_edges;
        const auto& mat = view.material_params;

        // init force and hessian
        Vec3 force = Vec3::Zero();
        Mat3 hessian = Mat3::Zero();

        // accumulate inertia force and hessian
        force += - (pos - inertia_pos) / (inv_mass * dt * dt);
        hessian += Mat3::Identity() / (inv_mass * dt * dt);

        // accumulate stvk triangle element force and hessian
        for (uint32_t f = face_adjacency.begin(vtex_id); f < face_adjacency.end(vtex_id); ++f ) {
            const auto pack = face_adjacency.incidents[f];
            const auto face_id = AdjacencyCSR::unpack_id(pack);
            const auto _order = AdjacencyCSR::unpack_order(pack);
            const auto& _face = view.tris[face_id];
            accumulate_stvk_triangle_force_hessian(view.pos, mat, _face, _order, force, hessian);
        }
        // accumulate bending element force and hessian
        for (uint32_t e = edge_adjacency.begin(vtex_id); e < edge_adjacency.end(vtex_id); ++e) {
            const auto pack = edge_adjacency.incidents[e];
            const auto edge_id = AdjacencyCSR::unpack_id(pack);
            const auto _order = AdjacencyCSR::unpack_order(pack);
            const auto& _edge = view.edges[edge_id];
            accumulate_dihedral_angle_based_bending_force_hessian(view.pos, mat, _edge, _order, force, hessian);
        }

        // solve linear system and update result
        const auto delta_x = TY::SolveSPDOrRegularize(hessian, force);
        pos += delta_x;
    }
}

void VBDSolver::update_velocity(SimView& view, const float dt) {
    for (size_t i = 0; i < view.pos.size(); ++i) {
        if (view.fixed[i]) {
            view.vel[i] = Vec3::Zero();
            continue;
        }
        view.vel[i] = (view.pos[i] - view.prev_pos[i]) / dt;
    }
}




























