//
// Created by 徐天焱 on 2025/11/5.
//

#ifndef TAIYI_RENDERHELPER_H
#define TAIYI_RENDERHELPER_H



#include <vector>
#include "raylib.h"
#include "Model.h" // MModel, State, MeshInfo, range, triangle
#include "Types.h" // Vec3, VertexID

class RenderHelper final {
public:
    RenderHelper() = default;
    ~RenderHelper() { Shutdown(); }

    RenderHelper(const RenderHelper&) = delete;
    RenderHelper& operator=(const RenderHelper&) = delete;

    void BindModel(const MModel& model);

    // 每帧：若拓扑版本变化 -> 自动重建；否则仅更新 VBO
    void Update(const State& state);

    void Draw() const;               // Draw 不需要 state
    void Shutdown();                 // 释放 GPU 资源
    [[nodiscard]] bool Ready() const { return ready_; }

private:
    struct RenderMesh {
        MeshInfo info{}; // 保存 range（值拷贝，避免 mesh_infos vector 重分配导致引用/指针风险）
        Model model{};   // raylib owning
        bool valid = false;
    };

private:
    void Rebuild();                  // 用 model_ 重建所有 GPU mesh/model
    void UpdateDynamic(const State& state);

    void CheckBound() const;

    static void FillPositionsXYZ(const State& state,
                                 size_t particle_begin,
                                 size_t particle_count,
                                 float* dst_xyz);

    static void ComputeNormalsXYZ(const MModel& model,
                                  const State& state,
                                  range tri_range,
                                  size_t particle_begin,
                                  size_t particle_count,
                                  float* dst_nxyz);

    static void BuildIndicesU16(const MModel& model,
                               range tri_range,
                               size_t particle_begin,
                               size_t particle_count,
                               unsigned short* dst_indices);

private:
    const MModel* model_ = nullptr;

    uint64_t built_topology_version_ = 0;
    bool ready_ = false;

    std::vector<RenderMesh> meshes_;
};




#endif //TAIYI_RENDERHELPER_H
