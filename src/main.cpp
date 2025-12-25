#include <functional>
#include <iostream>
#include <memory>
#include <vector>
#include <raylib.h>
#include <rlgl.h>

#include "World.h"
#include "RenderHelper.hpp"
#include "BroadPhase.h"
#include "MeshBuilder.h"
#include "Mesh.h"
#include "Types.h"
#include "CameraController.h"
#include "ShaderManager.h"


int main(){

    constexpr int screenWidth = 1920;
    constexpr int screenHeight = 1080;


    SetConfigFlags(FLAG_MSAA_4X_HINT); // 抗锯齿
    InitWindow(screenWidth, screenHeight, "tinyVBD");

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
            // ReframeOrbitToModels(orbitCam, models, 1.2f);
        }

        RefreshCameraTransform(orbitCam);

        // step simulation and update


        // ordinary rendering
        BeginDrawing();
        ClearBackground(RAYWHITE);
        DrawRectangleGradientV(
            0, 0, GetScreenWidth(), GetScreenHeight(),
            Color{ 45, 93, 125, 255 },   // 顶部：明显更亮的蓝灰
            Color{ 10, 12, 16, 255 }    // 底部：深色
        );


        BeginMode3D(orbitCam.camera);
        // 1. 画地板


        // 2. 画布料
        rlDisableBackfaceCulling();
        /*for (auto& m : models) {
            DrawModel(m, Vector3Zero(), 1.0f, DARKBLUE);
        }*/
        rlEnableBackfaceCulling();

        DrawAxisGizmo(0.8f);

        // DrawGrid(10, 1.0f); // 辅助网格可以留着，但它是线框，没阴影

        EndMode3D();
        EndDrawing();
    }

    // shaderManager.UnloadAll();
    /*for (auto& m : models) {
        UnloadModel(m);             // release gpu resource
    }*/
    // UnloadModel(floor);
    CloseWindow();

    return 0;
}
