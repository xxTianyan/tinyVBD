#include <vector>
#include <raylib.h>
#include <rlgl.h>
#include <rlImGui.h>
#include <imgui.h>

#include "RenderHelper.hpp"
#include "Mesh.h"
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
            // DrawModelWires(m, Vector3Zero(), 1.0f, DARKBLUE);
        }
        rlEnableBackfaceCulling();
        // draw coordinates
        // DrawAxisGizmo(0.8f);


        // adj test
        const auto& mesh = falling_cloth.m_world->meshes[0];
        const auto& tris = mesh->m_tris;
        const auto& adj_info = mesh->adjacencyInfo;
        const VertexId v = 177;
        for (const auto& t: tris) {
            auto v1 = ToRayVec(mesh->pos[t.vertices[0]]);
            auto v2 = ToRayVec(mesh->pos[t.vertices[1]]);
            auto v3 = ToRayVec(mesh->pos[t.vertices[2]]);
            DrawLine3D(v1, v2, MAROON);
            DrawLine3D(v2, v3, MAROON);
            DrawLine3D(v3, v1, MAROON);
        }
        DrawSphere(ToRayVec(mesh->pos[v]), 0.01, BLUE);
        for (uint32_t f = adj_info.vertex_faces.begin(v); f < adj_info.vertex_faces.end(v); ++f) {
            auto pack_id = adj_info.vertex_faces.incidents[f];
            auto face_id = AdjacencyCSR::unpack_id(pack_id);
            auto the_tri = mesh->m_tris[face_id];
            auto v1 = ToRayVec(mesh->pos[the_tri.vertices[0]]);
            auto v2 = ToRayVec(mesh->pos[the_tri.vertices[1]]);
            auto v3 = ToRayVec(mesh->pos[the_tri.vertices[2]]);
            DrawTriangle3D(v3, v2, v1, RAYWHITE);
        }


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
