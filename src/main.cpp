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


#define SHADOWMAP_RESOLUTION 2048

static RenderTexture2D LoadShadowmapRenderTexture(int width, int height) {
    RenderTexture2D target = { 0 };
    target.id = rlLoadFramebuffer();
    target.texture.width = width;
    target.texture.height = height;

    if (target.id > 0) {
        rlEnableFramebuffer(target.id);
        target.depth.id = rlLoadTextureDepth(width, height, false);
        target.depth.width = width;
        target.depth.height = height;
        target.depth.format = 19;
        target.depth.mipmaps = 1;
        rlFramebufferAttach(target.id, target.depth.id, RL_ATTACHMENT_DEPTH, RL_ATTACHMENT_TEXTURE2D, 0);
        if (rlFramebufferComplete(target.id)) TraceLog(LOG_INFO, "FBO: Shadowmap created successfully");
        rlDisableFramebuffer();
    }
    return target;
}

static void UnloadShadowmapRenderTexture(RenderTexture2D target) {
    if (target.id > 0) rlUnloadFramebuffer(target.id);
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
    constexpr size_t demo_id = 0;

    SetConfigFlags(FLAG_MSAA_4X_HINT); // 抗锯齿
    InitWindow(screenWidth, screenHeight, "tinyVBD");

    // Define the camera to look into our 3d world
    Camera3D camera = {0};
    camera.position = Vector3{ 1.5f, 1.5f, 1.5f }; // Camera position
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
    floor.materials[0].maps[MATERIAL_MAP_DIFFUSE].color = LIGHTGRAY;

    // shader
    Shader sh = LoadShader("../shaders/shadow_blinn.vs", "../shaders/shadow_blinn.fs");
    sh.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(sh, "viewPos");
    int lightDirLoc = GetShaderLocation(sh, "lightDir");
    int lightVPLoc  = GetShaderLocation(sh, "lightVP");
    int shadowMapLoc = GetShaderLocation(sh, "shadowMap");

    // set shader
    for (auto& m : models) m.materials[0].shader = sh;
    floor.materials[0].shader = sh;

    // shadow map
    RenderTexture2D shadowMap = LoadShadowmapRenderTexture(SHADOWMAP_RESOLUTION, SHADOWMAP_RESOLUTION);

    // light
    Vector3 lightPos = Vector3Normalize((Vector3){ 1.5f, 1.5f, 1.5f }); // 光照位置
    Camera3D lightCam = { 0 };
    lightCam.position = Vector3Scale(lightPos, 1.0f); // 光源位置
    lightCam.target = Vector3Zero();
    lightCam.projection = CAMERA_ORTHOGRAPHIC; // 定向光用正交投影
    lightCam.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    lightCam.fovy = 4.0f; // 正交视野宽度 (覆盖场景范围)


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

        // render light
        Matrix lightView;
        Matrix lightProj;
        Matrix lightViewProj;

        BeginTextureMode(shadowMap);
        ClearBackground(RAYWHITE); // 深度图初始化为最远
        BeginMode3D(lightCam);
        // 绘制场景中所有需要投射阴影的物体

        // 绘制布料 (注意：这一步不需要高级 Shader，只需要深度)
        // 这里的关键技巧：Pass 1 不需要绑定我们的 shadow shader
        // Raylib 的默认行为会写入深度缓冲，这对我们足够了
        // 如果你想更严谨，可以用一个空的 shader，但直接画通常也行
        for (auto& m : models) {
            DrawModel(m, Vector3Zero(), 1.0f, WHITE);
        }

        // 记录光照矩阵用于 Pass 2
        lightView = rlGetMatrixModelview();
        lightProj = rlGetMatrixProjection();
        EndMode3D();
        EndTextureMode();

        lightViewProj = MatrixMultiply(lightView, lightProj);

        // ordinary rendering
        BeginDrawing();
        ClearBackground(RAYWHITE);
        // update shader uniforms
        Vector3 camPos = camera.position;
        SetShaderValue(sh, sh.locs[SHADER_LOC_VECTOR_VIEW], &camPos, SHADER_UNIFORM_VEC3);
        SetShaderValue(sh, lightDirLoc, &lightPos, SHADER_UNIFORM_VEC3);
        SetShaderValueMatrix(sh, lightVPLoc, lightViewProj);

        // 激活 Shadow Map 纹理
        int slot = 10;
        rlActiveTextureSlot(slot);
        rlEnableTexture(shadowMap.depth.id);

        // 【关键修复】直接使用 shadowMapLoc 变量
        rlSetUniform(shadowMapLoc, &slot, SHADER_UNIFORM_INT, 1);

        BeginMode3D(camera);
        // 1. 画地板 (接收阴影)


        // 2. 画布料 (接收阴影 + 自投影)
        rlDisableBackfaceCulling();
        for (auto& m : models) {
            DrawModel(m, Vector3Zero(), 1.0f, ORANGE);
        }
        rlEnableBackfaceCulling();

        DrawModel(floor, Vector3Zero(), 1.0f, WHITE);

        // 3. 画光源位置示意
        DrawSphere(lightCam.position, 0.2f, YELLOW);
        // DrawGrid(10, 1.0f); // 辅助网格可以留着，但它是线框，没阴影

        EndMode3D();

        DrawTextureEx(shadowMap.depth, (Vector2){10, 10}, 0.0f, 0.1f, WHITE);
        // DrawText("Shadow Map Debug View", 10, 120, 10, GRAY);

        EndDrawing();
    }

    UnloadShader(sh);
    for (auto& m : models) {
        UnloadModel(m);             // release gpu resource
    }
    CloseWindow();

    return 0;
}