#include <iostream>
#include "Mesh.h"
#include "Types.h"
#include <functional>
#include <memory>

#include "World.h"
#include "RenderHelper.hpp"
#include "BroadPhase.h"


/* for initialize */
namespace {
    auto gravity = Vec3(0, 0, -9.81f);
    World world(gravity);
    // test if it can push to git on fysics computer

}

/************ demos ************/

static void demo1() {
    auto m = std::make_unique<mesh_on_cpu>();
    ParseMSH("../assets/bunny.msh", m.get());
    world.Add(std::move(m));
}

std::function<void()> demos[] = {demo1};

static void InitDemo(size_t index) {
    world.Clear();
    demos[index]();
}

int main(){

    constexpr int screenWidth = 1280;
    constexpr int screenHeight = 720;

    InitWindow(screenWidth, screenHeight, "tinyVBD");

    // Define the camera to look into our 3d world
    Camera3D camera = {0};
    camera.position = Vector3{ 1.5f, 1.5f, 1.5f }; // Camera position
    camera.target = Vector3{ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = Vector3{ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera projection type
    DisableCursor();                    // Limit cursor to relative movement inside the window

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

    InitDemo(0);
    const auto models = upload_all_models(world);

    // shader
    const Shader sh = LoadShader("../shaders/blinnphong_vs.glsl", "../shaders/blinnphong_fs.glsl");
    sh.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(sh, "viewPos");
    for (const auto& m : models) {
        m.materials[0].shader = sh;
    }
    // uniforms
    const int locLightPos     = GetShaderLocation(sh, "uLightPos");
    const int locLightColor   = GetShaderLocation(sh, "uLightColor");
    const int locAlbedoColor  = GetShaderLocation(sh, "uAlbedoColor");
    const int locAmbient      = GetShaderLocation(sh, "uAmbient");
    const int locShininess    = GetShaderLocation(sh, "uShininess");
    const int locSpecStrength = GetShaderLocation(sh, "uSpecStrength");
    // initial parameters
    constexpr Vector3 lightPos = { 0.0f, 3.0f, 0.0f };
    constexpr float   lightColor[3]   = { 1.0f, 1.0f, 1.0f };
    constexpr float   albedoColor[3]  = { 1.0f, 0.7f, 0.2f };   // model's adhere color
    constexpr float   ambient[3]      = { 0.35f, 0.35f, 0.35f };
    constexpr float   shininess       = 24.0f;
    constexpr float   specStrength    = 0.4f;
    // set shader
    SetShaderValue(sh, locLightPos,     &lightPos,     SHADER_UNIFORM_VEC3);
    SetShaderValue(sh, locLightColor,   lightColor,    SHADER_UNIFORM_VEC3);
    SetShaderValue(sh, locAlbedoColor,  albedoColor,   SHADER_UNIFORM_VEC3);
    SetShaderValue(sh, locAmbient,      ambient,       SHADER_UNIFORM_VEC3);
    SetShaderValue(sh, locShininess,    &shininess,    SHADER_UNIFORM_FLOAT);
    SetShaderValue(sh, locSpecStrength, &specStrength, SHADER_UNIFORM_FLOAT);

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
        world.Step(dt);
        UpdateModel(models, world.meshes);

        BeginDrawing();

        ClearBackground(RAYWHITE);
        BeginMode3D(camera);

        for (auto& m : models) {
            DrawModel(m, Vector3{0,0,0}, 1.0f, GRAY);
        }

        DrawGrid(10, 1.0);

        // DrawModelWires(model, Vector3{0,0,0}, 1.0f, BLACK);
        EndMode3D();

        EndDrawing();
    }
    for (auto& m : models) {
        UnloadModel(m);             // release gpu resource
    }
    CloseWindow();

    return 0;
}