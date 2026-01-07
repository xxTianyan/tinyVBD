//
// Created by xumiz on 2025/12/24.
//

#include "Sample.h"
#include "Application.h"
#include "VBDSolver.h"
void Sample::OnEnter(AppContext &ctx) {
    // set up scene
    CreateWorld(ctx);

    // set up floor
    CreateFloor(ctx);

    // upload mesh to gpu
    BuildRenderResources();

    // bind shaders
    BindShaders(ctx);

}

void Sample::OnExit([[maybe_unused]]AppContext &ctx) {
    // clean cpu resource
    CleanUp();

    // clean gpu resource, important!!
    DestroyRenderResources();
}

void Sample::Update([[maybe_unused]]AppContext &ctx) {
    if (ctx.paused) return;
    if (scene_ == nullptr) return;
    if (ctx.dt == 0)
        Step(1.0f / static_cast<float>(ctx.target_fps));
    else
        Step(ctx.dt);

    renderHelper_.Update(scene_->state_out_);
}

void Sample::Render([[maybe_unused]]AppContext &ctx) {
    BeginMode3D(ctx.orbitCam->camera);

    // floor
    if (RenderHelper::IsModelValid(floor_)) {
        DrawModel(floor_, Vector3{0,0,0}, 1.0f, WHITE);
    }

    // scene models
    renderHelper_.Draw();

    EndMode3D();
}

void Sample::DrawUI([[maybe_unused]]AppContext &ctx) {
    // no panel on default
}

void Sample::CleanUp() {
    solver_.reset();
    scene_.reset();
}

void Sample::BuildRenderResources() {
    renderHelper_.BindModel(scene_->model_);
    renderHelper_.Update(scene_->state_out_);      // need update once manually in case app is paused and pass sample update in main loop
}

void Sample:: DestroyRenderResources() {
    renderHelper_.Shutdown();
    RenderHelper::UnloadRLModelSafe(floor_);  // models that have no physical meanings is owned by sample itself
}

void Sample::CreateFloor([[maybe_unused]]AppContext& ctx) {

    if (ctx.shader_manager == nullptr)
        throw std::runtime_error("No shader manager found");

    ctx.shader_manager->LoadShaderProgram("floor", "../resources/shaders/floor.vs", "../resources/shaders/floor.fs");
    const auto floor_shader = ctx.shader_manager->Get("floor")->shader;
    ShaderManager::BindMatrices(floor_shader);
    ShaderManager::SetCommonShaderParams(floor_shader);
    const Mesh floor_mesh = GenMeshPlane(500.05f, 500.0f, 1, 1);
    floor_ = LoadModelFromMesh(floor_mesh);
    floor_.materials[0].shader = floor_shader;

    // locate uniform parameters
    const int tileScale = ShaderManager::CheckSetShaderLocation(floor_shader, "tileScale");
    const int lineWidth = ShaderManager::CheckSetShaderLocation(floor_shader, "lineWidth");
    const int baseAColor = ShaderManager::CheckSetShaderLocation(floor_shader, "baseAColor");
    const int baseBColor = ShaderManager::CheckSetShaderLocation(floor_shader, "baseBColor");
    const int lineColor = ShaderManager::CheckSetShaderLocation(floor_shader, "lineColor");
    const int roughness = ShaderManager::CheckSetShaderLocation(floor_shader, "roughness");
    const int bumpStrength = ShaderManager::CheckSetShaderLocation(floor_shader, "bumpStrength");
    const int fogDensity = ShaderManager::CheckSetShaderLocation(floor_shader, "fogDensity");
    const int fogColor = ShaderManager::CheckSetShaderLocation(floor_shader, "fogColor");


    // Floor appearance
    constexpr Vector3 floorFogColor  = { 0.10f, 0.13f, 0.17f }; // 用于 fogColor（线性，偏蓝灰）
    constexpr Vector3 floorBaseACol  = { 0.08f, 0.085f, 0.09f }; // 底色A：深
    constexpr Vector3 floorBaseBCol  = { 0.13f, 0.135f, 0.14f }; // 底色B：浅（噪声混合）
    constexpr Vector3 floorLineColor   = { 0.20f, 0.205f, 0.215f}; // 网格线颜色（别太亮）

    constexpr float floorRough  = 0.55f;                // 越大越哑光，反光越弱
    constexpr float floorBumpStr     = 0.22f;                // 微起伏：增强高级感（太大像橡皮泥）
    constexpr float floorFogDensity  = 0.015f;               // 雾：让地板“无穷远”+聚光更明显
    constexpr float floorTileScale   = 2.0f;                 // 网格密度：越大格子越小
    constexpr float floorLineWidth   = 0.035f;               // 网格线宽：越大线越明显

    SetShaderValue(floor_shader, tileScale, &floorTileScale, SHADER_UNIFORM_FLOAT);
    SetShaderValue(floor_shader, lineWidth, &floorLineWidth, SHADER_UNIFORM_FLOAT);
    SetShaderValue(floor_shader, baseAColor, &floorBaseACol, SHADER_UNIFORM_VEC3);
    SetShaderValue(floor_shader, baseBColor, &floorBaseBCol, SHADER_UNIFORM_VEC3);
    SetShaderValue(floor_shader, lineColor, &floorLineColor, SHADER_UNIFORM_VEC3);

    SetShaderValue(floor_shader, roughness, &floorRough, SHADER_UNIFORM_FLOAT);
    SetShaderValue(floor_shader, bumpStrength, &floorBumpStr, SHADER_UNIFORM_FLOAT);
    SetShaderValue(floor_shader, fogDensity, &floorFogDensity, SHADER_UNIFORM_FLOAT);
    SetShaderValue(floor_shader, fogColor, &floorFogColor, SHADER_UNIFORM_VEC3);
}


/*my_Sample::Sample() {
    m_world = std::make_unique<Scene>(Vec3{0.0f, -9.81f, 0.0f});
    m_shader_manager = std::make_unique<ShaderManager>();
    m_solver = std::make_unique<VBDSolver>(10);
}

void my_Sample::CleanUp() {
    auto UnloadModelIfLoaded = [](Model& m) {
        if ((m.meshCount > 0) && (m.meshes != nullptr)) {
            UnloadModel(m);
            m = Model{};
        }
    };
    UnloadModelIfLoaded(m_floor);
    for (auto& m : m_models) {
        UnloadModelIfLoaded(m);
    }
    UnloadModelIfLoaded(m_floor);
    m_shader_manager->UnloadAll();
}

void my_Sample::Step(const float dt) {
    m_world->InitStep();
    auto& meshes = m_world->meshes;
    for (size_t mesh_id = 0; mesh_id < meshes.size(); mesh_id++) {
        SimView view = m_world->MakeSimView(mesh_id);
        // make inertia step
        VBDSolver::forward_step(view, dt);
        // iter newton step
        for (size_t iter = 0; iter < 20; iter++)
            VBDSolver::solve(view, dt);
        VBDSolver::update_velocity(view, dt);
    }
}*/



