#include <vector>
#include <raylib.h>
#include <rlgl.h>
#include <rlImGui.h>
#include <imgui.h>

#include "RenderHelper.hpp"
#include "Mesh.h"
#include "CameraController.h"
#include "Sample.h"

int main(){

    constexpr int screenWidth = 800;  // 1920
    constexpr int screenHeight = 600; // 1080

    SetConfigFlags(FLAG_MSAA_4X_HINT); // anti alias
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    InitWindow(screenWidth, screenHeight, "tinyVBD");
    // initialize rlImGui
    rlImGuiSetup(true);
    // Monitor
    PerformanceMonitor perfMonitor;

    // set sample
    HangingCloth falling_cloth;
    falling_cloth.CreateWorld();
    falling_cloth.CreateFloor();
    falling_cloth.BindShaders();

    // get reference of necessary component
    auto& models = falling_cloth.m_models;
    const auto& shader_manager = falling_cloth.m_shader_manager;
    bool& isPaused = falling_cloth.isPaused;

    // camera
    OrbitCamera orbitCam = CreateOrbitCamera(Vector3{ 1.5f, 0.0f, 0.0f }, Vector3{ 0.0f, 0.0f, 0.0f });

    SetTargetFPS(60);                   // run at 60 frames-per-second


    // main loop
    while (!WindowShouldClose())        // Detect window close button or ESC key
    {
        const float dt = GetFrameTime() == 0? 1.0f/60 : GetFrameTime();
        perfMonitor.Update(dt);
        // not move camera when mouse is on imgui window
        if (!ImGui::GetIO().WantCaptureMouse) {
            auto [x, y] = GetMouseDelta();
            const float wheel = GetMouseWheelMove();
            UpdateOrbitCameraInput(orbitCam, Vector2{(float)x, (float)y}, wheel,
                                   IsMouseButtonDown(MOUSE_RIGHT_BUTTON), IsMouseButtonDown(MOUSE_MIDDLE_BUTTON));
        }

        if (IsKeyPressed(KEY_Z)) {
            ReframeOrbitToModels(orbitCam, models, 1.2f);
        }
        RefreshCameraTransform(orbitCam);

        if (IsKeyPressed(KEY_SPACE)) isPaused = !isPaused;

        // update shader
        const auto& viewPos = orbitCam.camera.position;
        shader_manager->UpdateViewPos("cloth", viewPos);
        shader_manager->UpdateViewPos("floor", viewPos);

        // step simulation and update model
        if (!isPaused) {
            // falling_cloth.Step(dt);
        }
        UpdateModel(models, falling_cloth.m_world->meshes);

        // ordinary rendering
        BeginDrawing();
        ClearBackground(RAYWHITE);
        DrawRectangleGradientV(
            0, 0, GetScreenWidth(), GetScreenHeight(),
            Color{ 45, 93, 125, 255 },   // 顶部：明显更亮的蓝灰
            Color{ 10, 12, 16, 255 }    // 底部：深色
        );

        BeginMode3D(orbitCam.camera);
        // draw floor
        DrawModel(falling_cloth.m_floor, Vector3Zero(), 1.0f, WHITE);
        // draw cloth
        rlDisableBackfaceCulling();
        for (const auto& m : models) {
            // DrawModel(m, Vector3Zero(), 1.0f, DARKBLUE);
            DrawModelWires(m, Vector3Zero(), 1.0f, DARKBLUE);
        }
        rlEnableBackfaceCulling();
        // draw coordinates
        DrawAxisGizmo(0.8f);
        // DrawGrid(10, 1.0f);
        EndMode3D();

        rlImGuiBegin();
        perfMonitor.Draw();
        rlImGuiEnd();

        EndDrawing();
    }
    falling_cloth.CleanUp();
    rlImGuiShutdown();
    CloseWindow();

    return 0;
}
