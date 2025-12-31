//
// Created by tianyan on 12/22/25.
//

#include "VBDDynamics.h"
#include <iostream>

static Vec3 SolveSPDOrRegularize(Mat3 H, const Vec3& f, const double eps = 1e-9) {
    using Scalar = typename Mat3::Scalar;

    // 强制对称化：用 Mat3 的标量类型，避免 double/float 推导问题
    constexpr auto half = static_cast<Scalar>(0.5);
    H = (half * (H + H.transpose())).eval();

    // 先尝试 LLT
    Eigen::LLT<Mat3> llt;
    llt.compute(H);
    if (llt.info() == Eigen::Success) {
        return llt.solve(f);
    }

    // 正则化：H + (eps*scale) * I   ——只平移对角，保持对称
    const Scalar diag_max = H.diagonal().cwiseAbs().maxCoeff();
    const Scalar scale    = std::max<Scalar>(static_cast<Scalar>(1), diag_max);

    H.diagonal().array() += static_cast<Scalar>(eps) * scale;

    // 再试一次 LLT
    llt.compute(H);
    if (llt.info() == Eigen::Success) {
        return llt.solve(f);
    }

    // 兜底：用 LDLT（对不完全 SPD 更鲁棒）
    Eigen::LDLT<Mat3> ldlt;
    ldlt.compute(H);
    return ldlt.solve(f);
}

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
    const Vec3 f0 = (x1 - x0) * DmInv00 + (x2 - x0) * DmInv01;
    const Vec3 f1 = (x1 - x0) * DmInv10 + (x2 - x0) * DmInv11;

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
        force = Vec3::Zero();
        H = Mat3::Zero();
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
    const Vec3 delta_force = -dpsi_dx;

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
   const Mat3 delta_hessian = df0_dx_sq * d2E_dF2_00 + df1_dx_sq * d2E_dF2_11 +
       df0_df1_cross * (d2E_dF2_01 + d2E_dF2_01.transpose());


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

    // make inertia step
    VBDSolver::forward_step(view, dt);

    // solve
    for (VertexID vtex_id = 0; vtex_id < num_nodes; ++vtex_id) {
        const auto& inv_mass = view.inv_mass[vtex_id];
        if (inv_mass == 0) continue;

        auto& pos = view.pos[vtex_id];
        auto& vel = view.vel[vtex_id];
        const auto& inertia_pos = view.inertia_pos[vtex_id];
        const auto& prev_pos = view.prev_pos[vtex_id];
        // const auto& edge_adjacency = view.adj.vertex_edges;
        const auto& face_adjacency = view.adj.vertex_faces;
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

        // solve linear system and update result
        const auto delta_x = SolveSPDOrRegularize(hessian, force);
        pos += delta_x;
        vel = (pos - prev_pos) / dt;
    }

}




























