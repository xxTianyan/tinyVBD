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

    static void accumulate_stvk_triangle_force_hessian(std::span<const Vec3> pos, const MMaterial& mat,
                                        const triangle& face, uint32_t vtex_order, Vec3& force, Mat3& H);

    static void accumulate_dihedral_angle_based_bending_force_hessian(std::span<const Vec3> pos, const MMaterial& mat,
                                        const edge& e, uint32_t vtex_order, Vec3& force, Mat3& H);

private:
    int num_iters;
    std::vector<Vec3> inertia;
    std::vector<Vec3> prev_pos;
    std::vector<ForceElementAdjacencyInfo> adjacencyInfo;
};

struct Node {
    VertexID idx = INVALID_VERTEX_ID;
};





#endif //TAIYI_VBDDYNAMICS_H
