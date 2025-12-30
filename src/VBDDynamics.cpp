//
// Created by tianyan on 12/22/25.
//

#include "VBDDynamics.h"

#include <iostream>

inline Vec3 SolveSPDOrRegularize(Mat3 H, const Vec3& f, const double eps = 1e-9) {
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


void VBDSolver::forward_step(SimView &view, const float dt) {
    const size_t num_nodes = view.pos.size();
    for (size_t i = 0; i < num_nodes; ++i) {
        view.prev_pos[i] = view.pos[i];  // record previous pos
        view.inertia_pos[i] = view.pos[i] + dt * view.vel[i] + dt * dt * view.accel[i];
    }
}

void VBDSolver::solve(SimView& view, const float dt) {
    const size_t num_nodes = view.pos.size();

    // make inertia step
    VBDSolver::forward_step(view, dt);

    // solve
    for (size_t i = 0; i < num_nodes; ++i) {
        const auto& inv_mass = view.inv_mass[i];
        if (inv_mass == 0) continue;

        auto& pos = view.pos[i];
        auto& vel = view.vel[i];
        const auto& inertia_pos = view.inertia_pos[i];
        const auto& prev_pos = view.prev_pos[i];

        // init force and hessian
        Vec3 force = Vec3::Zero();
        Mat3 hessian = Mat3::Zero();

        // accumulate inertia force and hessian
        force += - (pos - inertia_pos) / (inv_mass * dt * dt);
        hessian += I3() / (inv_mass * dt * dt);

        // accumulate potential force and hessian

        // solve linear system and update result
        const auto delta_x = SolveSPDOrRegularize(hessian, force);
        pos += delta_x;
        vel = (pos - prev_pos) / dt;
    }

}




























