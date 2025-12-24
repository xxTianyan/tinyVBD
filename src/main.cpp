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

struct DemoDefinition
{
    std::function<void(World &)> buildWorld;
    std::function<void(std::vector<Model> &, ShaderManager &)> bindShaders;
    std::function<void(ShaderManager &, const Vector3 &, const Vector3 &)> updateLighting;
};

/************ demos ************/

static DemoDefinition Demo0()
{
    DemoDefinition demo{};
    demo.buildWorld = [](World &world) {
        auto m = std::make_unique<mesh_on_cpu>();
        ParseMSH("../assets/bunny.msh", m.get());
        world.Add(std::move(m));
    };
    demo.bindShaders = [](std::vector<Model> &, ShaderManager &) {};
    demo.updateLighting = [](ShaderManager &, const Vector3 &, const Vector3 &) {};
    return demo;
}

static DemoDefinition Demo1()
{
    DemoDefinition demo{};
    demo.buildWorld = [](World &world) {
        auto m = std::make_unique<mesh_on_cpu>();
        MeshBuilder::BuildBox(m.get(), 1.0f, 1.0f, 1.0f);
        world.Add(std::move(m));
    };
    demo.bindShaders = [](std::vector<Model> &, ShaderManager &) {};
    demo.updateLighting = [](ShaderManager &, const Vector3 &, const Vector3 &) {};
    return demo;
}

static DemoDefinition Demo2()
{
    DemoDefinition demo{};
    demo.buildWorld = [](World &world) {
        auto m = std::make_unique<mesh_on_cpu>();
        MeshBuilder::BuildCloth(m.get(), 1.0f, 2.0f, 10, 20, Vec3{0.0f, 2.0f, 0.0f});
        world.Add(std::move(m));
    };
    demo.bindShaders = [](std::vector<Model> &models, ShaderManager &manager) {
        manager.LoadClothShader();
        manager.ApplyShaderToModels(models, "cloth");
    };
    demo.updateLighting = [](ShaderManager &manager, const Vector3 &lightDir, const Vector3 &viewPos) {
        manager.UpdateLighting("cloth", lightDir, viewPos);
    };
    return demo;
}

static DemoDefinition demos[] = {Demo0(), Demo1(), Demo2()};

static const DemoDefinition &InitDemo(size_t index)
{
    world.Clear();
    const DemoDefinition &demo = demos[index];
    demo.buildWorld(world);
    return demo;
}

int main(){

    constexpr int screenWidth = 1280;
    constexpr int screenHeight = 720;
    constexpr size_t demo_id = 2;

    SetConfigFlags(FLAG_MSAA_4X_HINT); // 抗锯齿
    InitWindow(screenWidth, screenHeight, "tinyVBD");

    OrbitCamera orbitCam = CreateOrbitCamera(Vector3{ 1.5f, 0.0f, 0.0f }, Vector3{ 0.0f, 0.0f, 0.0f });

    SetTargetFPS(60);                   // run at 60 frames-per-second

    const DemoDefinition &demo = InitDemo(demo_id);
    auto models = upload_all_models(world);

    // create floor model
    Model floor = LoadModelFromMesh(GenMeshPlane(10.0f, 10.0f, 10, 10));
    floor.materials[0].maps[MATERIAL_MAP_DIFFUSE].color = RAYWHITE;

    ShaderManager shaderManager;
    demo.bindShaders(models, shaderManager);

    Vector3 lightDir = Vector3Normalize(Vector3{-0.7f, -1.0f, -0.45f});


    // auto adjcentTets = BuildNodeTetAdj(world.meshes[0]->size(), world.meshes[0]->m_tets_local);


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

        // step simulation and update
        world.Step(dt);
        UpdateModel(models, world.meshes);

        // ordinary rendering
        BeginDrawing();
        ClearBackground(RAYWHITE);
        demo.updateLighting(shaderManager, lightDir, orbitCam.camera.position);

        BeginMode3D(orbitCam.camera);
        // 1. 画地板
        DrawModel(floor, Vector3Zero(), 1.0f, WHITE);

        // 2. 画布料
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

    shaderManager.UnloadAll();
    for (auto& m : models) {
        UnloadModel(m);             // release gpu resource
    }
    UnloadModel(floor);
    CloseWindow();

    return 0;
}
