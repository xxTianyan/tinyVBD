//
// Created by 徐天焱 on 2025/11/11.
//

#ifndef TAIYI_VBDDYNAMICS_H
#define TAIYI_VBDDYNAMICS_H

#include <span>
#include <Types.h>
#include "AdjacencyCSR.hpp"
#include "Scene.h"
#include "ISolver.h"
#include "MaterialParams.hpp"


struct DebugTriggerConfig {
    bool enabled = true;
    bool freeze_on_trigger = true;

    // 触发条件（可按你单位调）
    float dx_limit_scale = 0.05f;     // dx_limit = dx_limit_scale * avg_edge_length
    float J_min = 0.2f;               // min det(F) below this triggers
    bool trigger_on_first_contact = true; // 一旦出现任何穿透就触发
    bool trigger_on_inversion = true; // current signed volume <= 0 触发
    bool trigger_on_nan = true;       // 任意 NaN/Inf 触发
};

struct DebugFrameStats {
    // 全局统计
    float minJ = std::numeric_limits<float>::infinity();
    size_t minJ_tet = size_t(-1);

    float minSignedVol = std::numeric_limits<float>::infinity();
    size_t minVol_tet = size_t(-1);

    float minAbsVol = std::numeric_limits<float>::infinity();
    size_t minAbsVol_tet = size_t(-1);

    float minRestVol = std::numeric_limits<float>::infinity();
    size_t minRestVol_tet = size_t(-1);

    float maxPenetration = 0.0f;
    size_t maxPen_vtx = size_t(-1);

    // 触发点信息
    size_t trigger_vertex = size_t(-1);
    float trigger_dx_norm = 0.0f;
    float trigger_dx_limit = 0.0f;
    float trigger_pen = 0.0f;

    // 你也可以加 frame_id / iter_id
    size_t frame_id = 0;
};

class VBDSolver final : public ISolver {

public:
    explicit VBDSolver(const MModel* model, const int num_iters, const MMaterial& material = default_cloth())
        : model_(model), num_iters(num_iters), material_(material) {}
    ~VBDSolver() override = default;

    void Init() override;

    void Step(State& state_in, State& state_out, float dt) override;

    void forward_step(State& state_in, float dt);

    void solve_serial(State& state_in, State& state_out, float dt) const;

    void update_velocity(State& stat_out, float dt) const;

    static void accumulate_stvk_triangle_force_hessian(std::span<const Vec3> pos, const MMaterial& mat,
                                        const triangle& face, uint32_t vtex_order, Vec3& force, Mat3& H);

    static void accumulate_stvk_triangle_force_hessian_serial(std::span<const Vec3> pos, const MMaterial& mat,
                                    const triangle& face, uint32_t vtex_order, Vec3& force, Mat3& H);

    static void accumulate_dihedral_angle_based_bending_force_hessian(std::span<const Vec3> pos, const MMaterial& mat,
                                        const edge& e, uint32_t vtex_order, Vec3& force, Mat3& H);

    static void accumulate_dihedral_angle_based_bending_force_hessian_serial(std::span<const Vec3> pos, const MMaterial& mat,
                                    const edge& e, uint32_t vtex_order, Vec3& force, Mat3& H);

    static void accumulate_neo_hookean_tetrahedron_force_hessian(std::span<const Vec3> pos, const MMaterial& mat,
                                        const tetrahedron& tet, uint32_t vtex_order, Vec3& force, Mat3& H);

    // debug part
    bool DebugPaused() const noexcept { return debug_pause_; }
    void ClearDebugPause() noexcept { debug_pause_ = false; }
    void SetDebugTriggerConfig(const DebugTriggerConfig& cfg) { debug_cfg_ = cfg; }
    const DebugFrameStats& LastDebugStats() const noexcept { return last_debug_stats_; }

    DebugTriggerConfig debug_cfg_{};
    mutable bool debug_pause_ = false;
    mutable DebugFrameStats last_debug_stats_{};

    DebugFrameStats ComputeDebugStats(std::span<const Vec3> pos, float ground_y = 0.0f) const;
    void DumpDebugStats(const DebugFrameStats& s) const;

private:

    void BuildAdjacencyInfo();

    const MModel*  model_;

    int num_iters;

    MMaterial material_;

    std::vector<Vec3> inertia_;
    std::vector<Vec3> prev_pos_;

    ForceElementAdjacencyInfo adjacency_info_;

    // temporary
    std::vector<char> surface_vertices;

    uint64_t topology_version_ = 0;
};







#endif //TAIYI_VBDDYNAMICS_H
