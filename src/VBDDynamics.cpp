//
// Created by tianyan on 12/22/25.
//

#include "VBDDynamics.h"


void VBDSolver::forward_step(SimView &view, const float dt) {
    const size_t num_nodes = view.pos.size();
    for (size_t i = 0; i < num_nodes; ++i) {
        view.inertia_pos[i] = view.pos[i] + dt * view.vel[i] + dt * dt * view.accel[i];
    }
}

void VBDSolver::solve(SimView& view, const float dt) {
    const size_t num_nodes = view.pos.size();

    // make inertia step
    VBDSolver::forward_step(view, dt);

    // solve
    for (size_t i = 0; i < num_nodes; ++i) {
        view.prev_pos[i] = view.inertia_pos[i];  // for test
        view.vel[i] = (view.prev_pos[i] - view.pos[i]) / dt;
        view.pos[i] = view.prev_pos[i];
    }

}




























