#include <raylib.h>
#include <rlgl.h>
#include <rlImGui.h>
#include <imgui.h>

#include "RenderHelper.hpp"
#include "CameraController.h"
#include "Sample.h"

inline Vector3 ToRayVec(Vec3& v_pos) {
    return Vector3{v_pos.x(), v_pos.y(), v_pos.z()};
};

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

    // camera
    OrbitCamera orbitCam = CreateOrbitCamera(Vector3{ 1.5f, 0.0f, 0.0f }, Vector3{ 0.0f, 0.0f, 0.0f });

    SetTargetFPS(120);                   // run at 60 frames-per-second


    // main loop
    while (!WindowShouldClose())        // Detect window close button or ESC key
    {
        const float dt = GetFrameTime() == 0? 1.0f/120 : GetFrameTime();
        perfMonitor.Update(dt);
        // not move camera when mouse is on imgui window
        if (!ImGui::GetIO().WantCaptureMouse) {
            auto [x, y] = GetMouseDelta();
            const float wheel = GetMouseWheelMove();
            UpdateOrbitCameraMouse(orbitCam, Vector2{(float)x, (float)y}, wheel,
                                   IsMouseButtonDown(MOUSE_RIGHT_BUTTON), IsMouseButtonDown(MOUSE_MIDDLE_BUTTON));
        }

        if (!ImGui::GetIO().WantCaptureKeyboard) {
            UpdateOrbitCameraKeyboard(orbitCam, dt, 4.5f);
        }


        RefreshCameraTransform(orbitCam);



        // update shader
        const auto& viewPos = orbitCam.camera.position;




        // rendering
        BeginDrawing();
        ClearBackground(RAYWHITE);
        DrawRectangleGradientV(
            0, 0, GetScreenWidth(), GetScreenHeight(),
            Color{ 45, 93, 125, 255 },   // 顶部：明显更亮的蓝灰
            Color{ 10, 12, 16, 255 }    // 底部：深色
        );

        BeginMode3D(orbitCam.camera);
        // draw floor

        // draw cloth
        rlDisableBackfaceCulling();
        /*for (const auto& m : models) {
            DrawModel(m, Vector3Zero(), 1.0f, DARKBLUE);
            // DrawModelWires(m, Vector3Zero(), 1.0f, DARKBLUE);
        }*/
        rlEnableBackfaceCulling();
        // draw coordinates
        DrawAxisGizmo(0.8f);

        EndMode3D();

        rlImGuiBegin();
        perfMonitor.Draw();
        rlImGuiEnd();

        EndDrawing();
    }

    rlImGuiShutdown();
    CloseWindow();

    return 0;
}
