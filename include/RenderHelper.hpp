//
// Created by 徐天焱 on 2025/11/5.
//

#ifndef RENDERHELPER_H
#define RENDERHELPER_H

#include <iostream>
#include <limits>
#include <stdexcept>
#include <raylib.h>
#include <vector>
#include <imgui.h>

#include "Mesh.h"
#include "World.h"

static void DrawAxisGizmo(float length = 0.5f) {
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

// turn mesh_on_cpu into GPU Mesh + Model (dynamic)
static Model upload_model_from_cpu_mesh(mesh_on_cpu* M) {
    Mesh gmsh = {0};
    gmsh.vertexCount   = static_cast<int>(M->size());
    gmsh.triangleCount = static_cast<int>(M->m_surface_tris.size() / 3);

    const std::vector<float> vertices = assemble_vertices(M);
    // copy vertices
    gmsh.vertices = static_cast<float *>(MemAlloc(vertices.size() * sizeof(float)));
    std::memcpy(gmsh.vertices, vertices.data(), vertices.size() * sizeof(float));
    // copy indices
    gmsh.indices = static_cast<unsigned short *>(MemAlloc(M->m_surface_tris.size() * sizeof(unsigned short)));
    std::memcpy(gmsh.indices, M->m_surface_tris.data(), M->m_surface_tris.size() * sizeof(unsigned short));

    if (World::RayNormal) {
        const std::vector<float> normals  = ComputeNormal(M);
        gmsh.normals = static_cast<float *>(MemAlloc(normals.size() * sizeof(float)));
        std::memcpy(gmsh.normals, normals.data(), normals.size() * sizeof(float));
    }

    UploadMesh(&gmsh, /*dynamic=*/true);

    // raylib will hold the mesh
    const Model model = LoadModelFromMesh(gmsh);
    return model;
}

// put every mesh to a model vector
static std::vector<Model> upload_all_models(const World& world) {
    std::vector<Model> models;
    models.reserve(world.meshes.size());
    for (auto& up : world.meshes) {
        models.push_back(upload_model_from_cpu_mesh(up.get()));
    }
    return models;
}

// update gpu mesh according to cpu mesh
static void UpdateModel(const std::vector<Model> &models, const std::vector<MeshPtr>& Meshes) {
    if (models.size() != Meshes.size()) throw std::runtime_error("Model size mismatch");
    std::vector<float> vertices;
    std::vector<float> normals;
    for (size_t i = 0; i < models.size(); ++i) {
        const Mesh& old_mesh = models[i].meshes[0];
        vertices = assemble_vertices(Meshes[i].get());
        normals = ComputeNormal(Meshes[i].get());
        const size_t dataSize = vertices.size() * sizeof(float);
        if (dataSize > static_cast<size_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("Data too large for int parameter");
        }
        UpdateMeshBuffer(old_mesh, 0, vertices.data(), static_cast<int>(dataSize),0);
        UpdateMeshBuffer(old_mesh, 2, normals.data(), static_cast<int>(dataSize), 0);
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
        float frameTimeMs = dt * 1000.0f;
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
            ImGui::PlotLines("##FrameTimes", m_frameTimes.data(), HISTORY_SIZE, m_offset, nullptr, 0.0f, 33.3f, ImVec2(220, 50));
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

#endif //RENDERHELPER_H
