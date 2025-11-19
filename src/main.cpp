#include <iostream>
#include "Mesh.h"
#include "Types.h"
#include <functional>
#include "World.h"
#include "RenderHelper.hpp"
#include "BroadPhase.h"


/* for initialize */
namespace {
    auto gravity = Vec3(0, 0, -9.81f);
    World world(gravity);

}

/************  demos ************/

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
        radius > 0.001f ? radius : 3.0f, // 防止零半径
        0.0035f,   // 旋转灵敏度（鼠标像素 -> 弧度）
        0.12f,     // 滚轮每格的比例缩放
        0.0025f    // 中键平移灵敏度
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
    constexpr float   albedoColor[3]  = { 1.0f, 0.7f, 0.2f };   // 你的固有色
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


    // main loop
    while (!WindowShouldClose())        // Detect window close button or ESC key
    {

        UpdateCamera(&camera, CAMERA_FREE);
        // 在你的更新循环里加：
        if (IsKeyPressed(KEY_Z)) {
            ReframeToModel(&camera, &orbit, models, 1.2f); // Z：重置到屏幕中央并装满
        }

        // 鼠标增量
        auto [x, y] = GetMouseDelta();
        const float wheel = GetMouseWheelMove();

        // 右键拖动：绕 target 旋转
        if (IsMouseButtonDown(MOUSE_RIGHT_BUTTON)) {
            orbit.yaw   -= x * orbit.rotSens;
            orbit.pitch -= y * orbit.rotSens;
            ClampPitch(&orbit.pitch);
        }

        // 滚轮：缩放（改变半径）
        if (fabsf(wheel) > 0.0f) {
            const float factor = 1.0f - wheel * orbit.zoomSens; // 正轮缩小半径（拉近）
            // 限制缩放范围
            float newR = orbit.radius * factor;
            if (newR < 0.2f) newR = 0.2f;        // 最小距离
            if (newR > 100.0f) newR = 100.0f;    // 最大距离
            orbit.radius = newR;
        }

        // 中键拖动：平移 target（沿相机右/上方向）
        if (IsMouseButtonDown(MOUSE_MIDDLE_BUTTON)) {
            // 根据当前 yaw/pitch 求相机的前、右、上向量
            const Vector3 viewDir = Vector3Normalize(SphericalToCartesian(1.0f, orbit.yaw, orbit.pitch));
            const Vector3 right   = Vector3Normalize(Vector3CrossProduct(viewDir, camera.up));
            const Vector3 up      = Vector3Normalize(Vector3CrossProduct(right, viewDir));

            // 像素位移映射到世界位移：右移为 +x，中键上拖为 +y
            const Vector3 delta = Vector3Add(Vector3Scale(right, -x * orbit.panSens * orbit.radius),
                                       Vector3Scale(up,    y * orbit.panSens * orbit.radius));

            camera.target = Vector3Add(camera.target, delta);
        }

        // 根据球坐标回写 position
        camera.position = Vector3Add(camera.target,
                                     SphericalToCartesian(orbit.radius, orbit.yaw, orbit.pitch));

        // 维持世界上方向
        camera.up = (Vector3){0,1,0};

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