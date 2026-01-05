//
// Created by 徐天焱 on 2025/11/5.
//

#ifndef TAIYI_RENDERHELPER_H
#define TAIYI_RENDERHELPER_H

#include <iostream>
#include <limits>
#include <stdexcept>
#include <raylib.h>
#include <vector>
#include <imgui.h>

#include "Mesh.h"
#include "Scene.h"

static void DrawAxisGizmo(const float length = 0.5f) {
    constexpr Vector3 origin = { 0.0f, 0.0f, 0.0f };
    const Vector3 x = { length, 0.0f, 0.0f };
    const Vector3 y = { 0.0f, length, 0.0f };
    const Vector3 z = { 0.0f, 0.0f, length };

    DrawLine3D(origin, x, RED);
    DrawLine3D(origin, y, GREEN);
    DrawLine3D(origin, z, BLUE);

    const float radius = length * 0.03f;
    DrawSphereWires(x, radius, 1, 6, RED);
    DrawSphereWires(y, radius, 1, 6, GREEN);
    DrawSphereWires(z, radius, 1, 6, BLUE);
}

static void FillVerticesXYZ(const MeshState& state, float* dst_xyz) {
    const size_t n = state.pos.size();
    for (size_t i = 0; i < n; ++i) {
        const Vec3& p = state.pos[i];
        dst_xyz[i * 3 + 0] = p.x();
        dst_xyz[i * 3 + 1] = p.y();
        dst_xyz[i * 3 + 2] = p.z();
    }
}

// compute real-time vertex normal and upload to gpu later
static void FillNormalsXYZ(const MeshModel& model, const MeshState& state, float* dst_nxyz) {
    const size_t n = model.size();
    std::memset(dst_nxyz, 0, n * 3 * sizeof(float));

    const auto& idx = model.surface_tris;
    if (idx.size() % 3 != 0) {
        throw std::runtime_error("surface_tris size is not multiple of 3");
    }

    const size_t I = idx.size();
    for (size_t t = 0; t < I; t += 3) {
        const auto i0 = static_cast<size_t>(idx[t + 0]);
        const auto i1 = static_cast<size_t>(idx[t + 1]);
        const auto i2 = static_cast<size_t>(idx[t + 2]);

        if (i0 >= n || i1 >= n || i2 >= n) {
            throw std::runtime_error("surface_tris index out of range");
        }

        const Vec3& p0 = state.pos[i0];
        const Vec3& p1 = state.pos[i1];
        const Vec3& p2 = state.pos[i2];

        const Vec3 e1 = p1 - p0;
        const Vec3 e2 = p2 - p0;
        const Vec3 fn = e1.cross(e2); // 面法线（面积加权）

        dst_nxyz[i0 * 3 + 0] += fn.x();
        dst_nxyz[i0 * 3 + 1] += fn.y();
        dst_nxyz[i0 * 3 + 2] += fn.z();

        dst_nxyz[i1 * 3 + 0] += fn.x();
        dst_nxyz[i1 * 3 + 1] += fn.y();
        dst_nxyz[i1 * 3 + 2] += fn.z();

        dst_nxyz[i2 * 3 + 0] += fn.x();
        dst_nxyz[i2 * 3 + 1] += fn.y();
        dst_nxyz[i2 * 3 + 2] += fn.z();
    }

    // normalize per-vertex
    for (size_t i = 0; i < n; ++i) {
        const float nx = dst_nxyz[i * 3 + 0];
        const float ny = dst_nxyz[i * 3 + 1];
        const float nz = dst_nxyz[i * 3 + 2];

        if (const float len2 = nx * nx + ny * ny + nz * nz; len2 > 1e-30f) {
            const float inv = 1.0f / std::sqrt(len2);
            dst_nxyz[i * 3 + 0] = nx * inv;
            dst_nxyz[i * 3 + 1] = ny * inv;
            dst_nxyz[i * 3 + 2] = nz * inv;
        } else {
            dst_nxyz[i * 3 + 0] = 0.0f;
            dst_nxyz[i * 3 + 1] = 1.0f;
            dst_nxyz[i * 3 + 2] = 0.0f;
        }
    }
}

// upload all mesh info into GPU Mesh + Model (dynamic)
static Model upload_model_from_cpu_mesh(const MeshModel& model, const MeshState& state) {

    Mesh gmsh{};
    const size_t num_verts = state.pos.size();
    const size_t num_idx   = model.surface_tris.size();

    if (num_idx % 3 != 0) {
        throw std::runtime_error("surface_tris size is not multiple of 3");
    }

    gmsh.vertexCount   = to_int_checked(num_verts, "vertexCount");
    gmsh.triangleCount = to_int_checked(num_idx / 3, "triangleCount");

    // ---- vertices ----
    const size_t v_float_count = num_verts * 3;
    gmsh.vertices = static_cast<float*>(MemAlloc(v_float_count * sizeof(float)));
    if (!gmsh.vertices) throw std::runtime_error("MemAlloc failed: vertices");
    FillVerticesXYZ(state, gmsh.vertices);

    // ---- indices (raylib Mesh.indices is unsigned short*) ----
    if (num_verts > static_cast<size_t>(std::numeric_limits<unsigned short>::max()) + 1ull) {
        throw std::runtime_error("Too many vertices for 16-bit indices (raylib Mesh.indices is unsigned short).");
    }

    static_assert(sizeof(unsigned short) == sizeof(uint16_t),
              "raylib Mesh.indices is expected to be 16-bit (unsigned short).");

    gmsh.indices = static_cast<unsigned short*>(MemAlloc(num_idx * sizeof(unsigned short)));
    if (!gmsh.indices) throw std::runtime_error("MemAlloc failed: indices");

    // copy from model.surface_tris
    std::memcpy(gmsh.indices, model.surface_tris.data(), num_idx * sizeof(uint16_t));

    // ---- normals (optional) ----
    if (Scene::RayNormal) {
        gmsh.normals = static_cast<float*>(MemAlloc(v_float_count * sizeof(float)));
        if (!gmsh.normals) throw std::runtime_error("MemAlloc failed: normals");
        FillNormalsXYZ(model, state, gmsh.normals);
    }

    UploadMesh(&gmsh, /*dynamic=*/true);

    // Attention：LoadModelFromMesh will take gmsh memory（raylib convention），
    // do not release gmsh.vertices/indices/normals maunally
    const Model rl_model = LoadModelFromMesh(gmsh);
    return rl_model;
}

// put every mesh to a model vector
static std::vector<Model> upload_all_models(const Scene& scene) {
    const auto& models = scene.MeshModels();
    const auto& states = scene.MeshStates();
    std::vector<Model> out;
    out.reserve(models.size());
    for (size_t i = 0; i < models.size();++i)
        out.push_back(upload_model_from_cpu_mesh(models[i], states[i]));
    return out;
}

// update gpu mesh according to cpu mesh
static void UpdateModel(const std::vector<Model>& gpu_models,
                        const std::vector<MeshModel>& cpu_models,
                        const std::vector<MeshState>& cpu_states) {
    if (gpu_models.size() != cpu_models.size() || gpu_models.size() != cpu_states.size()) {
        throw std::runtime_error("Model size mismatch");
    }

    std::vector<float> vertices; // 复用缓冲，避免每帧 malloc/free
    std::vector<float> normals;

    for (size_t i = 0; i < gpu_models.size(); ++i) {
        const MeshModel& mm = cpu_models[i];
        const MeshState& ms = cpu_states[i];

        const size_t n = mm.size();
        const size_t float_count = n * 3;

        // --- vertices ---
        vertices.resize(float_count);
        FillVerticesXYZ(ms, vertices.data());

        const size_t v_bytes = vertices.size() * sizeof(float);
        if (v_bytes > static_cast<size_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("Vertex buffer too large for int parameter");
        }

        // raylib: meshes[0] 是主 mesh
        const Mesh& gpu_mesh = gpu_models[i].meshes[0];

        // slot 0 = vertices
        UpdateMeshBuffer(gpu_mesh, 0, vertices.data(), static_cast<int>(v_bytes), 0);

        // --- normals (optional) ---
        if (Scene::RayNormal) {
            normals.resize(float_count);
            FillNormalsXYZ(mm, ms, normals.data());

            const size_t n_bytes = normals.size() * sizeof(float);
            if (n_bytes > static_cast<size_t>(std::numeric_limits<int>::max())) {
                throw std::runtime_error("Normal buffer too large for int parameter");
            }

            // slot 2 = normals
            UpdateMeshBuffer(gpu_mesh, 2, normals.data(), static_cast<int>(n_bytes), 0);
        }
    }
}

// imgui part

class PerformanceMonitor {
public:
    PerformanceMonitor() {
        m_frameTimes.resize(HISTORY_SIZE, 0.0f);
    }

    // 在主循环开始处调用，传入上一帧的时间 (GetFrameTime())
    void Update(const float dt) {
        // 1. 存入历史数据 (用于画波形图)
        // 将 dt 转换为毫秒
        const float frameTimeMs = dt * 1000.0f;
        m_frameTimes[m_offset] = frameTimeMs;
        m_offset = (m_offset + 1) % HISTORY_SIZE;

        // 2. 只有当计时器到达间隔时才更新显示的文字数字 (避免数字跳动过快看不清)
        m_timer += dt;
        if (m_timer >= UPDATE_INTERVAL) {
            m_timer = 0.0f;

            // 计算当前瞬时值
            m_displayFps = (dt > 0.0f) ? 1.0f / dt : 0.0f;
            m_displayMs = frameTimeMs;

            // 计算统计值 (遍历历史记录)
            float sum = 0.0f;
            float minT = 9999.0f;
            float maxT = 0.0f;
            int validCount = 0;

            for (float t : m_frameTimes) {
                if (t <= 0.0f) continue;
                sum += t;
                if (t < minT) minT = t;
                if (t > maxT) maxT = t;
                validCount++;
            }

            if (validCount > 0) {
                // FrameTime(ms) 越低，FPS 越高
                m_avgFps = 1000.0f / (sum / validCount);
                m_minFps = (maxT > 0.0f) ? 1000.0f / maxT : 0.0f; // 最慢的一帧决定了最低FPS
                m_maxFps = (minT > 0.0f) ? 1000.0f / minT : 0.0f; // 最快的一帧决定了最高FPS
            }
        }
    }

    // 在 rlImGuiBegin() 和 rlImGuiEnd() 之间调用
    void Draw() const {
        // --- 样式配置 ---
        // 设置半透明深色背景，看起来像专业的 Overlay
        ImGui::SetNextWindowBgAlpha(0.15f);

        // 窗口标志：无标题栏、自动大小、不可聚焦(不抢键盘输入)、禁止导航
        ImGuiWindowFlags windowFlags = ImGuiWindowFlags_NoDecoration |
                                       ImGuiWindowFlags_AlwaysAutoResize |
                                       ImGuiWindowFlags_NoFocusOnAppearing |
                                       ImGuiWindowFlags_NoNav;

        // 固定在左上角 (Padding 10px)
        constexpr float PAD = 0.0f;
        ImGui::SetNextWindowPos(ImVec2(PAD, PAD), ImGuiCond_Always);

        if (ImGui::Begin("PerfMonitor", nullptr, windowFlags)) {

            // 1. 顶部状态栏：渲染后端信息 (这里写死，如果有API可以获取动态的)
            ImGui::TextDisabled("Raylib 5.5 | V-Sync: ON");
            ImGui::Separator();

            // 2. 主 FPS 显示 (大字体)
            ImGui::SetWindowFontScale(1.8f);

            // 颜色动态变化：低帧率红色，高帧率绿色
            ImVec4 fpsColor = (m_displayFps < 30.0f) ? ImVec4(1.0f, 0.4f, 0.4f, 1.0f) : ImVec4(0.4f, 1.0f, 0.4f, 1.0f);
            float plotWidth = ImGui::GetContentRegionAvail().x;
            if (plotWidth < 200.0f) plotWidth = 220.0f;

            ImGui::TextColored(fpsColor, "%.0f FPS", m_displayFps);
            ImGui::SetWindowFontScale(1.0f); // 还原字体大小

            ImGui::SameLine();
            ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "(%.2f ms)", m_displayMs);

            ImGui::Spacing();

            // 3. 统计数据表格 (AVG / 1% LOW / MAX)
            if (ImGui::BeginTable("PerfStats", 3, ImGuiTableFlags_BordersInnerV)) {
                // 表头
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::TextDisabled("AVG");
                ImGui::TableSetColumnIndex(1); ImGui::TextDisabled("MIN");
                ImGui::TableSetColumnIndex(2); ImGui::TextDisabled("MAX");

                // 数据
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("%.0f", m_avgFps);
                ImGui::TableSetColumnIndex(1); ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.5f, 1.0f), "%.0f", m_minFps);
                ImGui::TableSetColumnIndex(2); ImGui::TextColored(ImVec4(0.5f, 1.0f, 1.0f, 1.0f), "%.0f", m_maxFps);

                ImGui::EndTable();
            }

            ImGui::Separator();

            // 4. 实时波形图 (Frame Times)
            ImGui::Spacing();
            ImGui::TextDisabled("Frame Time Graph (16ms target)");

            // PlotLines 参数解释：label, values, count, offset, overlay_text, scale_min, scale_max, graph_size
            // scale_max 设为 33.3ms (即 30FPS)，超过这个数值的波峰会被切断，适合观察流畅度
            ImGui::PlotLines(
                "##FrameTimes", m_frameTimes.data(), HISTORY_SIZE, m_offset, nullptr, 0.0f, 33.3f, ImVec2(plotWidth, 50));
        }
        ImGui::End();
    }

private:
    static constexpr int HISTORY_SIZE = 120;     // 记录最近120帧
    const float UPDATE_INTERVAL = 0.25f;     // 文字每0.25秒刷新一次

    std::vector<float> m_frameTimes;
    int m_offset = 0;
    float m_timer = 0.0f;

    // 缓存的显示数据
    float m_displayFps = 0.0f;
    float m_displayMs = 0.0f;
    float m_avgFps = 0.0f;
    float m_minFps = 0.0f;
    float m_maxFps = 0.0f;
};

#endif //TAIYI_RENDERHELPER_H
