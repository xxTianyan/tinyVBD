#include <raylib.h>
#include <iostream>
#include "Mesh.h"
#include "raymath.h"
#include "Types.h"


typedef struct OrbitCtrl {
    float yaw;        // 水平角 (弧度)
    float pitch;      // 俯仰角 (弧度)
    float radius;     // 目标点到相机的距离
    float rotSens;    // 旋转灵敏度
    float zoomSens;   // 缩放灵敏度
    float panSens;    // 平移灵敏度(像素->世界)
} OrbitCtrl;

static inline void ClampPitch(float *pitch) {
    constexpr float lim = DEG2RAD*89.0f;                 // 避免万向节锁
    if (*pitch >  lim) *pitch =  lim;
    if (*pitch < -lim) *pitch = -lim;
}

static inline Vector3 SphericalToCartesian(const float r, const float yaw, const float pitch) {
    const float cp = cosf(pitch);
    const float sp = sinf(pitch);
    const float cy = cosf(yaw);
    const float sy = sinf(yaw);
    return (Vector3){ r*cp*cy, r*sp, r*cp*sy };
}

// 用包围盒居中并设定合理半径
static inline void ReframeToModel(Camera3D *cam, OrbitCtrl *orbit, const Model &model, float margin)
{
    auto [min, max] = GetModelBoundingBox(model);
    const Vector3 center = {
        (min.x + max.x)*0.5f,
        (min.y + max.y)*0.5f,
        (min.z + max.z)*0.5f
    };
    // 用包围盒对角线的一半当“半径”，确保整盒子能装进视口
    const float radiusBox = 0.5f * Vector3Length(Vector3Subtract(max, min));
    const float fitDist   = radiusBox / tanf(DEG2RAD*cam->fovy*0.5f);

    orbit->radius = fitDist * (margin > 1.0f ? margin : 1.15f); // 留边 15%+
    cam->target   = center;

    // 按当前 yaw/pitch 重新放置相机
    const Vector3 offset = SphericalToCartesian(orbit->radius, orbit->yaw, orbit->pitch);
    cam->position  = Vector3Add(cam->target, offset);
    cam->up        = (Vector3){0,1,0};
}

int main(){
    constexpr int screenWidth = 1280;
    constexpr int screenHeight = 720;

    mesh_on_cpu cmesh;
    const std::vector<tetrahedron> tets = ParseMSH("../assets/bunny.msh", cmesh);
    const NodeTetAdj adj = buildNodeTetAdj(cmesh.size(), tets);
    const IndexBuffer tri_indices = BuildSurfaceTriangles(tets);

    const int vertexCount   = static_cast<int>(cmesh.size());
    const int triangleCount = static_cast<int>(tri_indices.size() / 3);

    Mesh gmesh = {0};
    gmesh.vertexCount = vertexCount;
    gmesh.triangleCount = triangleCount;

    InitWindow(screenWidth, screenHeight, "tinyVBD");

    // Define the camera to look into our 3d world
    Camera3D camera = {0};
    camera.position = (Vector3){ 1.5f, 1.5f, 1.5f }; // Camera position
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera projection type
    DisableCursor();                    // Limit cursor to relative movement inside the window

    // 计算初始球坐标（基于当前 position/target）
    Vector3 toCam = Vector3Subtract(camera.position, camera.target);
    float radius  = Vector3Length(toCam);
    float yaw     = atan2f(toCam.z, toCam.x);              // [-PI, PI]
    float pitch   = asinf(toCam.y / radius);               // [-PI/2, PI/2]
    ClampPitch(&pitch);

    OrbitCtrl orbit = {
        .yaw = yaw,
        .pitch = pitch,
        .radius = radius > 0.001f ? radius : 3.0f, // 防止零半径
        .rotSens  = 0.0035f,   // 旋转灵敏度（鼠标像素 -> 弧度）
        .zoomSens = 0.12f,     // 滚轮每格的比例缩放
        .panSens  = 0.0025f    // 中键平移灵敏度
    };

    SetTargetFPS(60);                   // run at 60 frames-per-second

    // create gmesh buffer
    const std::vector<float> normals = ComputeNormal(cmesh, tri_indices);
    const std::vector<float> vertices = assemble_vertices(cmesh);

    // copy vertices
    gmesh.vertices = static_cast<float *>(MemAlloc(vertices.size() * sizeof(float)));
    std::memcpy(gmesh.vertices, vertices.data(), vertices.size() * sizeof(float));
    // copy normal
    gmesh.normals = static_cast<float *>(MemAlloc(normals.size() * sizeof(float)));
    std::memcpy(gmesh.normals, normals.data(), normals.size() * sizeof(float));
    // copy indices
    gmesh.indices = static_cast<unsigned short *>(MemAlloc(tri_indices.size() * sizeof(unsigned short)));
    std::memcpy(gmesh.indices, tri_indices.data(), tri_indices.size() * sizeof(unsigned short));

    UploadMesh(&gmesh, /*dynamic=*/true);
    Model model = LoadModelFromMesh(gmesh);

    // shader
    Shader sh = LoadShader("../shaders/blinnphong_vs.glsl", "../shaders/blinnphong_fs.glsl");
    sh.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(sh, "viewPos");
    model.materials[0].shader = sh;
    // uniforms
    int locLightPos     = GetShaderLocation(sh, "uLightPos");
    int locLightColor   = GetShaderLocation(sh, "uLightColor");
    int locAlbedoColor  = GetShaderLocation(sh, "uAlbedoColor");
    int locAmbient      = GetShaderLocation(sh, "uAmbient");
    int locShininess    = GetShaderLocation(sh, "uShininess");
    int locSpecStrength = GetShaderLocation(sh, "uSpecStrength");
    // initial parameters
    Vector3 lightPos = { 0.0f, 3.0f, 0.0f };
    float   lightColor[3]   = { 1.0f, 1.0f, 1.0f };
    float   albedoColor[3]  = { 1.0f, 0.7f, 0.2f };   // 你的固有色
    float   ambient[3]      = { 0.35f, 0.35f, 0.35f };
    float   shininess       = 24.0f;
    float   specStrength    = 0.4f;
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
            ReframeToModel(&camera, &orbit, model, 1.2f); // Z：重置到屏幕中央并装满
        }

        // 鼠标增量
        auto [x, y] = GetMouseDelta();
        float wheel = GetMouseWheelMove();

        // 右键拖动：绕 target 旋转
        if (IsMouseButtonDown(MOUSE_RIGHT_BUTTON)) {
            orbit.yaw   -= x * orbit.rotSens;
            orbit.pitch -= y * orbit.rotSens;
            ClampPitch(&orbit.pitch);
        }

        // 滚轮：缩放（改变半径）
        if (fabsf(wheel) > 0.0f) {
            float factor = 1.0f - wheel * orbit.zoomSens; // 正轮缩小半径（拉近）
            // 限制缩放范围
            float newR = orbit.radius * factor;
            if (newR < 0.2f) newR = 0.2f;        // 最小距离
            if (newR > 100.0f) newR = 100.0f;    // 最大距离
            orbit.radius = newR;
        }

        // 中键拖动：平移 target（沿相机右/上方向）
        if (IsMouseButtonDown(MOUSE_MIDDLE_BUTTON)) {
            // 根据当前 yaw/pitch 求相机的前、右、上向量
            Vector3 viewDir = Vector3Normalize(SphericalToCartesian(1.0f, orbit.yaw, orbit.pitch));
            Vector3 right   = Vector3Normalize(Vector3CrossProduct(viewDir, camera.up));
            Vector3 up      = Vector3Normalize(Vector3CrossProduct(right, viewDir));

            // 像素位移映射到世界位移：右移为 +x，中键上拖为 +y
            Vector3 delta = Vector3Add(Vector3Scale(right, -x * orbit.panSens * orbit.radius),
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

        DrawModel(model, Vector3{0,0,0}, 1.0f, GRAY);
        DrawGrid(10, 1.0);

        // DrawModelWires(model, Vector3{0,0,0}, 1.0f, BLACK);
        EndMode3D();

        EndDrawing();
    }
    UnloadModel(model);               // release gpu resource
    CloseWindow();

    return 0;
}