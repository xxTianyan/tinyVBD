//
// Created by xumiz on 2026/1/1.
//

#include <stdexcept>
#include "Application.h"
#include "raylib.h"
#include "imgui.h"
#include "rlImGui.h"
#include "CameraController.h"


Application::Application(Desc desc) : desc_(desc) {
    // bind cxt
    ctx_.shader_manager = &shader_manager_;
    ctx_.orbitCam  = CreateOrbitCamera(Vector3{ 1.5f, 0.0f, 0.0f }, Vector3{ 0.0f, 0.0f, 0.0f });
}

Application::~Application() {
    if (IsWindowReady()) {
        ExitCurrentSample_();
        ShutdownWindowAndUI_();
    }
}

void Application::InitWindowAndUI_() const {
    // Window flags
    if (desc_.resizable) SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    InitWindow(desc_.width, desc_.height, desc_.title);
    SetTargetFPS(desc_.target_fps);

    // rlImGui init
    rlImGuiSetup(true); // true: enable docking if supported by your build
}

void Application::ShutdownWindowAndUI_() {
    rlImGuiShutdown();
    CloseWindow();
}

void Application::EnterSample_(const SampleId id) {
    current_id_ = id;
    current_ = registry_.Create(id);
    if (!current_) throw std::runtime_error("Create sample returned null");

    current_->OnEnter(ctx_);
}

void Application::ExitCurrentSample_() {
    if (current_) {
        current_->OnExit(ctx_);
        current_.reset();
    }
}

void Application::RequestSwitchSample_(const SampleId id) {
    pending_.type = PendingActionType::SwitchSample;
    pending_.target = id;
}

void Application::RequestReloadSample_() {
    pending_.type = PendingActionType::ReloadSample;
    pending_.target = current_id_;
}

void Application::PollHotkeys_() {
    /*// Reset
    if (IsKeyPressed(KEY_R)) {
        RequestReloadSample_();
    }
    // Pause
    if (IsKeyPressed(KEY_SPACE)) {
        pending_.type = PendingActionType::TogglePause;
    }
    // Reload shaders
    if (IsKeyPressed(KEY_F5)) {
        pending_.type = PendingActionType::ReloadShaders;
    }

    // 快捷切 sample（示例：1/2/3）
    if (IsKeyPressed(KEY_ONE))  RequestSwitchSample_(SampleId::ClothDrop);
    if (IsKeyPressed(KEY_TWO))  RequestSwitchSample_(SampleId::TetStVK);*/
}

void Application::DrawAppUI_() {

    ImGui::Begin("Application");

    ImGui::Text("Current: %s", current_ ? current_->Name() : "(none)");
    ImGui::Separator();

    // Pause / Reset / Reload shaders
    if (ImGui::Checkbox("Paused", &ctx_.paused)) {
    }
    if (ImGui::Button("Reset (Reload Sample) [R]")) {
        RequestReloadSample_();
    }
    if (ImGui::Button("Reload Shaders [F5]")) {
        pending_.type = PendingActionType::ReloadShaders;
    }

    ImGui::Separator();
    ImGui::Text("Switch Sample:");

    // pull down to choose sample
    static int current_index = 0;
    if (const auto& infos = registry_.Infos(); !infos.empty()) {
        // 同步 index（当你按热键切换后，UI 也要跟上）
        for (int i = 0; i < static_cast<int>(infos.size()); ++i) {
            if (infos[i].id == current_id_) { current_index = i; break; }
        }

        if (ImGui::BeginCombo("Sample", infos[current_index].display_name)) {
            for (int i = 0; i < static_cast<int>(infos.size()); ++i) {
                const bool selected = (i == current_index);
                if (ImGui::Selectable(infos[i].display_name, selected)) {
                    current_index = i;
                    RequestSwitchSample_(infos[i].id);
                }
                if (selected) ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }
    } else {
        ImGui::Text("No samples registered.");
    }

    ImGui::End();

    // Sample 自己的 UI
    if (current_) {
        current_->DrawUI(ctx_);
    }
}

void Application::ExecutePending_() {
    switch (pending_.type) {
        case PendingActionType::None:
            return;

        case PendingActionType::TogglePause:
            ctx_.paused = !ctx_.paused;
            break;

        case PendingActionType::ReloadShaders:
            // 注意：如果 Sample 内部缓存了 shader location/uniform handle，
            // 你需要让 Sample 在下一帧重新抓取，或把该逻辑放到 ShaderManager 封装里。
            // shader_manager_.ReloadAll();
            break;

        case PendingActionType::ReloadSample: {
            // 关键：销毁并重建同一个 sample（最稳 reset）
            ExitCurrentSample_();
            EnterSample_(pending_.target);
        } break;

        case PendingActionType::SwitchSample: {
            ExitCurrentSample_();
            EnterSample_(pending_.target);
        } break;
    }

    pending_.type = PendingActionType::None;
}

void Application::Run(const SampleId start_sample) {
    InitWindowAndUI_();

    // init sample
    EnterSample_(start_sample);

    while (!WindowShouldClose()) {
        ctx_.dt = GetFrameTime();

        PollHotkeys_();

        // Update
        if (current_ && !ctx_.paused) {
            current_->Update(ctx_);
        }

        // Render
        BeginDrawing();
        ClearBackground(RAYWHITE);

        if (current_) {
            current_->Render(ctx_);
        }

        // ImGui
        rlImGuiBegin();
        DrawAppUI_();
        rlImGuiEnd();

        EndDrawing();

        // delay excuting
        ExecutePending_();
    }

    ExitCurrentSample_();
    ShutdownWindowAndUI_();
}
