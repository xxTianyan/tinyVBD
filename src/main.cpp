#include <iostream>
#include <functional>
#include <memory>
#include <raylib.h>
#include <rlgl.h>

#include "World.h"
#include "RenderHelper.hpp"
#include "BroadPhase.h"
#include "MeshBuilder.h"
#include "Mesh.h"
#include "Types.h"

static void DrawAxisGizmo(float length = 0.5f) {
    const Vector3 origin = { 0.0f, 0.0f, 0.0f };
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


/* for initialize */
namespace {
    auto gravity = Vec3(0, 0, -9.81f);
    World world(gravity);
    // test if it can push to git on fysics computer

}

/************ demos ************/

static void demo0() {
    auto m = std::make_unique<mesh_on_cpu>();
    ParseMSH("../assets/bunny.msh", m.get());
    world.Add(std::move(m));
}

static void demo1() {
    auto m = std::make_unique<mesh_on_cpu>();
    MeshBuilder::BuildBox(m.get(), 1.0f ,1.0f, 1.0f);
    world.Add(std::move(m));
}

static void demo2() {
    auto m = std::make_unique<mesh_on_cpu>();
    MeshBuilder::BuildCloth(m.get(), 1.0f, 2.0f, 10, 20, Vec3{0.0f,2.0f,0.0f});
    world.Add(std::move(m));
}

std::function<void()> demos[] = {demo0, demo1, demo2};

static void InitDemo(size_t index) {
    world.Clear();
    demos[index]();
}

int main(){

    constexpr int screenWidth = 1280;
    constexpr int screenHeight = 720;
    constexpr size_t demo_id = 2;

    SetConfigFlags(FLAG_MSAA_4X_HINT); // 抗锯齿
    InitWindow(screenWidth, screenHeight, "tinyVBD");

    // Define the camera to look into our 3d world
    Camera3D camera = {0};
    camera.position = Vector3{ 1.5f, 0.0f, 0.0f }; // Camera position
    camera.target = Vector3{ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = Vector3{ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera projection type
    // DisableCursor();                    // Limit cursor to relative movement inside the window

    // 计算初始球坐标（基于当前 position/target）
    const Vector3 toCam = Vector3Subtract(camera.position, camera.target);
    const float radius  = Vector3Length(toCam);
    const float yaw     = atan2f(toCam.z, toCam.x);              // [-PI, PI]
    float pitch   = asinf(toCam.y / radius);               // [-PI/2, PI/2]
    ClampPitch(&pitch);

    OrbitCtrl orbit = {
        yaw,
        pitch,
        radius > 0.001f ? radius : 3.0f, // Prevent zero radius
        0.0035f,   // Rotation sensitivity
        0.12f,     // Wheel zoom step scale
        0.0025f    // Middle button pan sensitivity
    };

    SetTargetFPS(60);                   // run at 60 frames-per-second

    InitDemo(demo_id);
    auto models = upload_all_models(world);

    // create floor model
    Model floor = LoadModelFromMesh(GenMeshPlane(10.0f, 10.0f, 10, 10));
    floor.materials[0].maps[MATERIAL_MAP_DIFFUSE].color = RAYWHITE;

    // shader
    Shader sh = LoadShader("../shaders/shadow_blinn.vs", "../shaders/shadow_blinn.fs");
    sh.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(sh, "viewPos");
    int lightDirLoc = GetShaderLocation(sh, "lightDir");
    int lightVPLoc  = GetShaderLocation(sh, "lightVP");
    int shadowMapLoc = GetShaderLocation(sh, "shadowMap");

    // set shader
    for (auto& m : models) m.materials[0].shader = sh;

    // light
    Vector3 lightPos = { 6.0f, 6.0f, 6.0f }; // 将光源放高并远离，正交体能覆盖地面
    Vector3 lightDir = Vector3Normalize(Vector3Negate(lightPos));
    /*Camera3D lightCam = { 0 };
    lightCam.position = lightPos; // 光源位置
    lightCam.target = Vector3Zero();
    lightCam.projection = CAMERA_ORTHOGRAPHIC; // 定向光用正交投影
    lightCam.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    lightCam.fovy = 10.0f; // 正交视野宽度 (覆盖场景范围)*/


    // auto adjcentTets = BuildNodeTetAdj(world.meshes[0]->size(), world.meshes[0]->m_tets_local);


    // main loop
    while (!WindowShouldClose())        // Detect window close button or ESC key
    {
        float dt = GetFrameTime() == 0? 1.0f/60 : GetFrameTime();

        UpdateCamera(&camera, CAMERA_FREE);

        if (IsKeyPressed(KEY_Z)) {
            ReframeToModel(&camera, &orbit, models, 1.2f); // Z：reset view to center
        }

        auto [x, y] = GetMouseDelta();
        const float wheel = GetMouseWheelMove();

        if (IsMouseButtonDown(MOUSE_RIGHT_BUTTON)) {
            orbit.yaw   -= x * orbit.rotSens;
            orbit.pitch -= y * orbit.rotSens;
            ClampPitch(&orbit.pitch);
        }

        if (fabsf(wheel) > 0.0f) {
            const float factor = 1.0f - wheel * orbit.zoomSens;
            float newR = orbit.radius * factor;
            if (newR < 0.2f) newR = 0.2f;
            if (newR > 100.0f) newR = 100.0f;
            orbit.radius = newR;
        }

        if (IsMouseButtonDown(MOUSE_MIDDLE_BUTTON)) {
            const Vector3 viewDir = Vector3Normalize(SphericalToCartesian(1.0f, orbit.yaw, orbit.pitch));
            const Vector3 right   = Vector3Normalize(Vector3CrossProduct(viewDir, camera.up));
            const Vector3 up      = Vector3Normalize(Vector3CrossProduct(right, viewDir));

            const Vector3 delta = Vector3Add(Vector3Scale(right, -x * orbit.panSens * orbit.radius),
                                       Vector3Scale(up,    y * orbit.panSens * orbit.radius));

            camera.target = Vector3Add(camera.target, delta);
        }

        camera.position = Vector3Add(camera.target,
                                     SphericalToCartesian(orbit.radius, orbit.yaw, orbit.pitch));

        // keep with world up direction
        camera.up = {0,1,0};

        // step simulation and update
        world.Step(dt);
        UpdateModel(models, world.meshes);

        // ordinary rendering
        BeginDrawing();
        ClearBackground(RAYWHITE);
        // update shader uniforms
        Vector3 camPos = camera.position;
        SetShaderValue(sh, sh.locs[SHADER_LOC_VECTOR_VIEW], &camPos, SHADER_UNIFORM_VEC3);
        SetShaderValue(sh, lightDirLoc, &lightDir, SHADER_UNIFORM_VEC3);
        
        BeginMode3D(camera);
        // 1. 画地板 (接收阴影)
        DrawModel(floor, Vector3Zero(), 1.0f, WHITE);

        // 2. 画布料 (接收阴影 + 自投影)
        rlDisableBackfaceCulling();
        for (auto& m : models) {
            DrawModel(m, Vector3Zero(), 1.0f, ORANGE);
        }
        rlEnableBackfaceCulling();

        DrawAxisGizmo(0.8f);

        // DrawGrid(10, 1.0f); // 辅助网格可以留着，但它是线框，没阴影

        EndMode3D();
        EndDrawing();
    }

    UnloadShader(sh);
    for (auto& m : models) {
        UnloadModel(m);             // release gpu resource
    }
    CloseWindow();

    return 0;
}