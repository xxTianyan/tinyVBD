//
// Created by tianyan on 1/16/26.
//

#ifndef TAIYI_DBGUI_HPP
#define TAIYI_DBGUI_HPP

//
// Created by tianyan on 1/15/26.
// SolverDebuggerUI.hpp
// 依赖: imgui.h, TaiyiDebugger.hpp
//

#ifndef TAIYI_DEBUGGER_UI_HPP
#define TAIYI_DEBUGGER_UI_HPP

#include <string>
#include <ctime>
#include <sstream>
#include <iomanip>
#include "imgui.h"
#include "Debugger.hpp"

class SolverDebuggerUI {
public:
    // --- 主渲染函数 ---
    // 在你的 ImGui 循环中调用此函数
    // p_open: 可选的关闭窗口布尔指针
    void Render(SolverDebugger& debugger, bool* p_open = nullptr) {
        if (!ImGui::Begin("Solver Debugger Panel", p_open, ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::End();
            return;
        }

        DrawStatusAndControls(debugger);
        DrawManualDump(debugger);
        DrawTriggerInfo(debugger);
        DrawLiveStats(debugger);
        DrawSettings(debugger);

        ImGui::End();
    }

private:
    // --- 内部 UI 组件 ---

    void DrawStatusAndControls(SolverDebugger& debugger) {
        const auto& stats = debugger.get_current_stats();
        bool frozen = debugger.is_frozen();

        // 1. 状态灯
        ImGui::Text("Sim State: ");
        ImGui::SameLine();
        if (frozen) {
            ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "[ FROZEN ]");
            if (stats.trigger_reason != DebugErrorType::None) {
                ImGui::SameLine();
                ImGui::TextDisabled("(Triggered by Error)");
            }
        } else {
            ImGui::TextColored(ImVec4(0.4f, 1.0f, 0.4f, 1.0f), "[ RUNNING ]");
        }

        ImGui::Separator();

        // 2. VCR 控制器
        // Resume
        ImGui::BeginDisabled(!frozen);
        if (ImGui::Button("Resume (F5)")) {
            debugger.ui_continue();
        }
        ImGui::EndDisabled();

        ImGui::SameLine();

        // Pause
        ImGui::BeginDisabled(frozen);
        if (ImGui::Button(" Frozen ")) {
            debugger.ui_pause();
        }
        ImGui::EndDisabled();

        ImGui::SameLine();

        // Step
        if (ImGui::Button(" Step One ")) {
            debugger.ui_step_one();
        }
    }

    void DrawManualDump(SolverDebugger& debugger) {
        ImGui::Separator();
        ImGui::AlignTextToFramePadding();
        ImGui::Text("Snapshots:");
        ImGui::SameLine();

        if (ImGui::Button("Dump JSON History")) {
            std::string fname = GenerateTimestampFilename("debug_dump", ".json");
            debugger.dump_history_json(fname);
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Save current history ring buffer to JSON");
        }
    }

    void DrawTriggerInfo(SolverDebugger& debugger) {
        const auto& stats = debugger.get_current_stats();
        bool frozen = debugger.is_frozen();

        if (frozen && stats.trigger_reason != DebugErrorType::None) {
            ImGui::Separator();
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.6f, 0.6f, 1.0f));
            ImGui::Text("[!] CRASH REPORT");
            ImGui::PopStyleColor();

            if (ImGui::BeginTable("TriggerTable", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("Reason");
                ImGui::TableSetColumnIndex(1); ImGui::Text("%s", stats.trigger_msg.c_str());

                if (stats.trigger_element_id != static_cast<size_t>(-1)) {
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("Element ID");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%zu", stats.trigger_element_id);
                }

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("Value");
                ImGui::TableSetColumnIndex(1); ImGui::Text("%.6f", stats.trigger_val_1);

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("Force");
                ImGui::TableSetColumnIndex(1); ImGui::Text("(%.2f, %.2f, %.2f)",
                    stats.trigger_force.x(), stats.trigger_force.y(), stats.trigger_force.z());

                ImGui::EndTable();
            }
        }
    }

    void DrawLiveStats(const SolverDebugger& debugger) {
        if (!ImGui::CollapsingHeader("Live Statistics", ImGuiTreeNodeFlags_DefaultOpen)) return;

        const auto& stats = debugger.get_current_stats();

        ImGui::Text("Frame: %zu", stats.frame_id);

        // 辅助 lambda：显示带颜色的数值
        auto StatRow = [](const char* name, float val, float limit, bool is_min_limit) {
            bool bad = is_min_limit ? (val < limit) : (val > limit);
            ImGui::Text("%s:", name);
            ImGui::SameLine(280);
            if (bad) ImGui::TextColored(ImVec4(1, 0.3f, 0.3f, 1), "%.4f (BAD)", val);
            else     ImGui::Text("%.4f", val);
        };

        StatRow("Min Jacobian", stats.minJ, 0.1f, true);
        StatRow("Min Volume", stats.minSignedVol, 0.0f, true);
        StatRow("Max Penetration", stats.maxPenetration, 1e-3f, false);
        StatRow("Max Dx", stats.maxDx, 1.0f, false);
    }

    void DrawSettings(SolverDebugger& debugger) {
        if (!ImGui::CollapsingHeader("Debugger Config", ImGuiTreeNodeFlags_DefaultOpen)) return;

        DebugTriggerConfig& cfg = debugger.cfg_;

        ImGui::TextDisabled("Triggers (Auto-Freeze)");
        ImGui::Checkbox("NaN Detected", &cfg.break_on_nan);
        ImGui::Checkbox("Inverted (Vol<0)", &cfg.break_on_inversion);

        ImGui::Separator();
        ImGui::TextDisabled("Thresholds");

        // 使用 SliderFloat 而不是 DragFloat，可以防止用户拖出离谱的值
        ImGui::SliderFloat("Min Jacobian", &cfg.J_min_trigger, 0.0f, 0.5f);
        ImGui::SliderFloat("Max Penetration", &cfg.max_pen_trigger, 0.0f, 1.0f);
        ImGui::SliderFloat("Dx Limit Scale", &cfg.dx_limit_scale, 0.01f, 2.0f);
    }

    // --- 工具函数 ---
    std::string GenerateTimestampFilename(const std::string& prefix, const std::string& ext) {
        auto now = std::time(nullptr);
        auto tm = *std::localtime(&now);
        std::ostringstream oss;
        oss << prefix << "_" << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S") << ext;
        return oss.str();
    }
};

#endif // TAIYI_DEBUGGER_UI_HPP

#endif //TAIYI_DBGUI_HPP