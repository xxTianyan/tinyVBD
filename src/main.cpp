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
    Vector3 lightDir;
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
        manager.LoadFloorShader();
        manager.ApplyShaderToModels(models, "cloth");

        auto clothShader = manager.Get("cloth")->shader;
        clothShader.locs[SHADER_LOC_MATRIX_MVP]    = GetShaderLocation(clothShader, "mvp");
        clothShader.locs[SHADER_LOC_MATRIX_MODEL]  = GetShaderLocation(clothShader, "matModel");
        clothShader.locs[SHADER_LOC_MATRIX_NORMAL] = GetShaderLocation(clothShader, "matNormal");

        // 3) Register diffuse map sampler (texture0)
        clothShader.locs[SHADER_LOC_MAP_DIFFUSE]   = GetShaderLocation(clothShader, "texture0");

        // 4) Get custom uniform locations
        int loc_viewPos           = GetShaderLocation(clothShader, "viewPos");

        int loc_lightDir          = GetShaderLocation(clothShader, "lightDir");
        int loc_lightColor        = GetShaderLocation(clothShader, "lightColor");

        int loc_skyAmbientColor   = GetShaderLocation(clothShader, "skyAmbientColor");
        int loc_groundAmbientColor= GetShaderLocation(clothShader, "groundAmbientColor");
        int loc_ambientStrength   = GetShaderLocation(clothShader, "ambientStrength");

        int loc_shininess         = GetShaderLocation(clothShader, "shininess");
        int loc_specStrength      = GetShaderLocation(clothShader, "specStrength");
        int loc_wrapDiffuse       = GetShaderLocation(clothShader, "wrapDiffuse");

        int loc_exposure          = GetShaderLocation(clothShader, "exposure");

        // 5) Create material and assign shader
        Material clothMat = LoadMaterialDefault();
        clothMat.shader = clothShader;

        // Base color goes into colDiffuse automatically (raylib uses material.maps[ALBEDO].color)
        clothMat.maps[MATERIAL_MAP_ALBEDO].color = (Color){ 240, 200, 120, 255 };

        // If you have an albedo texture, set it here (optional)
        // clothMat.maps[MATERIAL_MAP_ALBEDO].texture = yourTexture;

        // 6) Set initial lighting/material parameters (dark scene preset)

        // Light direction computed from (lightPos -> target) each frame OR set once if static.
        // Here we set initial values; update in the frame loop too.
        Vector3 lightPos    = { 6.0f, 10.0f, 6.0f };
        Vector3 lightTarget = { 0.0f,  0.5f, 0.0f };
        lightDir    = Vector3Normalize(Vector3Subtract(lightTarget, lightPos));

        // HDR light intensity allowed (> 1.0)
        Vector3 lightColor  = { 1.0f, 1.0f, 1.0f };

        // Hemi ambient
        Vector3 skyAmbientColor    = { 0.10f, 0.12f, 0.18f };
        Vector3 groundAmbientColor = { 0.06f, 0.06f, 0.06f };
        float ambientStrength      = 1.0f;

        // Cloth-ish material
        float shininess    = 64.0f;
        float specStrength = 0.30f;
        float wrapDiffuse  = 0.25f;

        // Tonemap
        float exposure     = 1.6f;

        // 7) Upload uniforms (static-ish ones can be set once)
        SetShaderValue(clothShader, loc_lightDir,           &lightDir,           SHADER_UNIFORM_VEC3);
        SetShaderValue(clothShader, loc_lightColor,         &lightColor,         SHADER_UNIFORM_VEC3);

        SetShaderValue(clothShader, loc_skyAmbientColor,    &skyAmbientColor,    SHADER_UNIFORM_VEC3);
        SetShaderValue(clothShader, loc_groundAmbientColor, &groundAmbientColor, SHADER_UNIFORM_VEC3);
        SetShaderValue(clothShader, loc_ambientStrength,    &ambientStrength,    SHADER_UNIFORM_FLOAT);

        SetShaderValue(clothShader, loc_shininess,          &shininess,          SHADER_UNIFORM_FLOAT);
        SetShaderValue(clothShader, loc_specStrength,       &specStrength,       SHADER_UNIFORM_FLOAT);
        SetShaderValue(clothShader, loc_wrapDiffuse,        &wrapDiffuse,        SHADER_UNIFORM_FLOAT);

        SetShaderValue(clothShader, loc_exposure,           &exposure,           SHADER_UNIFORM_FLOAT);
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
    Model floor = LoadModelFromMesh(GenMeshPlane(500.0f, 500.0f, 1, 1));

    ShaderManager shaderManager;
    demo.bindShaders(models, shaderManager);

    // floor shader
    auto floorShader = shaderManager.Get("floor")->shader;

    floorShader.locs[SHADER_LOC_MATRIX_MVP]    = GetShaderLocation(floorShader, "mvp");
    floorShader.locs[SHADER_LOC_MATRIX_MODEL]  = GetShaderLocation(floorShader, "matModel");
    floorShader.locs[SHADER_LOC_MATRIX_NORMAL] = GetShaderLocation(floorShader, "matNormal");
    floor.materials[0].shader = floorShader;

    int loc_viewPos = GetShaderLocation(floorShader, "viewPos");
    int loc_lightDir = GetShaderLocation(floorShader, "lightDir");
    int loc_lightColor = GetShaderLocation(floorShader, "lightColor");
    int loc_sky = GetShaderLocation(floorShader, "skyAmbientColor");
    int loc_gnd = GetShaderLocation(floorShader, "groundAmbientColor");
    int loc_amb = GetShaderLocation(floorShader, "ambientStrength");
    int loc_time = GetShaderLocation(floorShader, "iTime");

    int loc_tile = GetShaderLocation(floorShader, "tileScale");
    int loc_lineW = GetShaderLocation(floorShader, "lineWidth");
    int loc_rough = GetShaderLocation(floorShader, "roughness");
    int loc_bump = GetShaderLocation(floorShader, "bumpStrength");
    int loc_fog = GetShaderLocation(floorShader, "fogDensity");
    int loc_exp = GetShaderLocation(floorShader, "exposure");

    Vector3 lightPos    = { 6.0f, 10.0f, 6.0f };
    Vector3 lightTarget = { 0.0f,  0.5f, 0.0f };
    lightDir    = Vector3Normalize(Vector3Subtract(lightTarget, lightPos));
    // HDR light intensity allowed (> 1.0)
    Vector3 lightColor  = { 1.0f, 1.0f, 1.0f };

    Vector3 skyAmb = { 0.10f, 0.12f, 0.18f };
    Vector3 gndAmb = { 0.05f, 0.05f, 0.055f };
    float ambientStrength = 0.06f;

    float tileScale = 5.0f;
    float lineWidth = 0.035f;
    float roughness = 0.50f;
    float bumpStrength = 0.22f;
    float fogDensity = 0.03f;
    float exposure = 1.6f;

    // set static-ish uniforms
    SetShaderValue(floorShader, loc_lightDir, &lightDir, SHADER_UNIFORM_VEC3);
    SetShaderValue(floorShader, loc_lightColor, &lightColor, SHADER_UNIFORM_VEC3);
    SetShaderValue(floorShader, loc_sky, &skyAmb, SHADER_UNIFORM_VEC3);
    SetShaderValue(floorShader, loc_gnd, &gndAmb, SHADER_UNIFORM_VEC3);
    SetShaderValue(floorShader, loc_amb, &ambientStrength, SHADER_UNIFORM_FLOAT);

    SetShaderValue(floorShader, loc_tile, &tileScale, SHADER_UNIFORM_FLOAT);
    SetShaderValue(floorShader, loc_lineW, &lineWidth, SHADER_UNIFORM_FLOAT);
    SetShaderValue(floorShader, loc_rough, &roughness, SHADER_UNIFORM_FLOAT);
    SetShaderValue(floorShader, loc_bump, &bumpStrength, SHADER_UNIFORM_FLOAT);
    SetShaderValue(floorShader, loc_fog, &fogDensity, SHADER_UNIFORM_FLOAT);
    SetShaderValue(floorShader, loc_exp, &exposure, SHADER_UNIFORM_FLOAT);

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
        DrawRectangleGradientV(
            0, 0, GetScreenWidth(), GetScreenHeight(),
            (Color){ 45, 93, 125, 255 },   // 顶部：明显更亮的蓝灰
            (Color){ 10, 12, 16, 255 }    // 底部：深色
        );
        demo.updateLighting(shaderManager, lightDir, orbitCam.camera.position);

        BeginMode3D(orbitCam.camera);
        // 1. 画地板
        DrawModel(floor, Vector3Zero(), 1.0f, WHITE);

        // 2. 画布料
        rlDisableBackfaceCulling();
        for (auto& m : models) {
            DrawModel(m, Vector3Zero(), 1.0f, DARKBLUE);
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
