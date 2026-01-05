//
// Created by 徐天焱 on 2025/11/11.
//

#ifndef TAIYI_VBDDYNAMICS_H
#define TAIYI_VBDDYNAMICS_H

#include <Types.h>
#include "Scene.h"

class VBDSolver {

public:
    explicit VBDSolver(const int num_iters) : num_iters(num_iters) {}
    ~VBDSolver() = default;

    static void forward_step(SimView& view, float dt);

    static void solve(SimView& view, float dt);

    static void update_velocity(SimView& view, float dt);

private:
    int num_iters;

};

struct Node {
    VertexID idx = INVALID_VERTEX_ID;
};





#endif //TAIYI_VBDDYNAMICS_H
