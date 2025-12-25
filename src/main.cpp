#include <functional>
#include <memory>
#include <vector>
#include <raylib.h>
#include <rlgl.h>


#include "RenderHelper.hpp"
#include "Mesh.h"
#include "CameraController.h"
#include "Sample.h"


int main(){

    constexpr int screenWidth = 1920;
    constexpr int screenHeight = 1080;

    SetConfigFlags(FLAG_MSAA_4X_HINT); // anti alias
    InitWindow(screenWidth, screenHeight, "tinyVBD");

    // set sample
    HangingCloth falling_cloth;
    falling_cloth.CreateWorld();
    falling_cloth.CreateFloor();
    falling_cloth.BindShaders();

    // get reference of necessary component
    auto& models = falling_cloth.m_models;
    auto& world = falling_cloth.m_world;
    const auto& shader_manager = falling_cloth.m_shader_manager;

    // camera
    OrbitCamera orbitCam = CreateOrbitCamera(Vector3{ 1.5f, 0.0f, 0.0f }, Vector3{ 0.0f, 0.0f, 0.0f });

    SetTargetFPS(60);                   // run at 60 frames-per-second


    // main loop
    while (!WindowShouldClose())        // Detect window close button or ESC key
    {
        float dt = GetFrameTime() == 0? 1.0f/60 : GetFrameTime();

        auto [x, y] = GetMouseDelta();
        const float wheel = GetMouseWheelMove();
        UpdateOrbitCameraInput(orbitCam, Vector2{static_cast<float>(x), static_cast<float>(y)}, wheel,
                               IsMouseButtonDown(MOUSE_RIGHT_BUTTON), IsMouseButtonDown(MOUSE_MIDDLE_BUTTON));

        if (IsKeyPressed(KEY_Z)) {
            ReframeOrbitToModels(orbitCam, models, 1.2f);
        }
        RefreshCameraTransform(orbitCam);

        const auto& viewPos = orbitCam.camera.position;

        // step simulation and update
        falling_cloth.m_world->Step(dt);
        // update shader
        shader_manager->UpdateViewPos("cloth", viewPos);
        shader_manager->UpdateViewPos("floor", viewPos);

        // ordinary rendering
        BeginDrawing();
        ClearBackground(RAYWHITE);
        DrawRectangleGradientV(
            0, 0, GetScreenWidth(), GetScreenHeight(),
            Color{ 45, 93, 125, 255 },   // 顶部：明显更亮的蓝灰
            Color{ 10, 12, 16, 255 }    // 底部：深色
        );
        DrawFPS(0, 0);

        BeginMode3D(orbitCam.camera);
        // draw floor
        DrawModel(falling_cloth.m_floor, Vector3Zero(), 1.0f, WHITE);

        // draw cloth
        rlDisableBackfaceCulling();
        for (const auto& m : models) {
            DrawModel(m, Vector3Zero(), 1.0f, DARKBLUE);
        }
        rlEnableBackfaceCulling();

        DrawAxisGizmo(0.8f);

        // DrawGrid(10, 1.0f);

        EndMode3D();
        EndDrawing();
    }

    CloseWindow();

    return 0;
}
