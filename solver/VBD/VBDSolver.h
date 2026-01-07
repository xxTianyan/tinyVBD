//
// Created by 徐天焱 on 2025/11/11.
//

#ifndef TAIYI_VBDDYNAMICS_H
#define TAIYI_VBDDYNAMICS_H

#include <span>
#include <Types.h>
#include "Scene.h"
#include "ISolver.h"
#include "MaterialParams.hpp"



class VBDSolver : public ISolver {

public:
    explicit VBDSolver(const int num_iters, const MMaterial& material = default_cloth())
        : num_iters(num_iters), material_(material) {}
    ~VBDSolver() = default;

    void Init(const Scene& scene) override;

    void Step(Scene& scene, float dt) override;

    void forward_step(Scene& scene, float dt);

    void solve(Scene& scene, float dt);

    void update_velocity(Scene& scene, float dt);

    static void accumulate_stvk_triangle_force_hessian(std::span<const Vec3> pos, const MMaterial& mat,
                                        const triangle& face, uint32_t vtex_order, Vec3& force, Mat3& H);

    static void accumulate_dihedral_angle_based_bending_force_hessian(std::span<const Vec3> pos, const MMaterial& mat,
                                        const edge& e, uint32_t vtex_order, Vec3& force, Mat3& H);

private:
    int num_iters;
    MMaterial material_;
    std::vector<Vec3> inertia_;
    std::vector<Vec3> prev_pos_;
    ForceElementAdjacencyInfo adjacency_info_;
    uint64_t topology_version_ = 0;
};







#endif //TAIYI_VBDDYNAMICS_H
