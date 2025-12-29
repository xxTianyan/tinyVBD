//
// Created by xumiz on 2025/12/24.
//

#include "Sample.h"
#include "MeshBuilder.h"
#include "RenderHelper.hpp"

Sample::Sample() {
    m_world = std::make_unique<World>(Vec3{0.0f, -9.81f, 0.0f});
    m_shader_manager = std::make_unique<ShaderManager>();
    m_solver = std::make_unique<VBDSolver>(10);
}

void Sample::CleanUp() {
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

void Sample::Step(const float dt) {
    m_world->InitStep();
    auto& meshes = m_world->meshes;
    for (auto& m:meshes) {
        SimView view = World::MakeSimView(*m);
        VBDSolver::solve(view, dt);
    }
}

void Sample::CreateFloor() {
    m_shader_manager->LoadShaderProgram("floor", "../shaders/floor.vs", "../shaders/floor.fs");
    const auto floor_shader = m_shader_manager->Get("floor")->shader;
    ShaderManager::BindMatrices(floor_shader);
    ShaderManager::SetCommonShaderParams(floor_shader);
    const Mesh floor_mesh = GenMeshPlane(500.05f, 500.0f, 1, 1);
    m_floor = LoadModelFromMesh(floor_mesh);
    m_floor.materials[0].shader = floor_shader;

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

void HangingCloth::CreateWorld() {
    m_world->ChangeGravity(Vec3{0.0f, -1.5f, 0.0f});
    auto m = std::make_unique<mesh_on_cpu>();
    MeshBuilder::BuildCloth(m.get(), 1.0f, 2.0f, 10, 20, Vec3{0.0f, 2.0f, 0.0f}, ClothOrientation::Horizontal);
    if (!m_world) throw std::runtime_error("m_world is empty pointer");
    m_world->Add(std::move(m));
    m_models = upload_all_models(*m_world);
}

void HangingCloth::BindShaders() const {
    m_shader_manager->LoadShaderProgram("cloth", "../shaders/cloth.vs", "../shaders/cloth.fs");
    const auto cloth_shader = m_shader_manager->Get("cloth")->shader;
    ShaderManager::BindMatrices(cloth_shader);
    ShaderManager::SetCommonShaderParams(cloth_shader);
    cloth_shader.locs[SHADER_LOC_MAP_DIFFUSE] = GetShaderLocation(cloth_shader, "texture0");
    auto m_cloth = m_models[0];
    m_cloth.materials[0].shader = cloth_shader;
    m_cloth.materials[0].maps[MATERIAL_MAP_ALBEDO].color = Color{230, 200, 160, 255};
    // m_cloth.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = clothAlbedoTex;          // 如果有贴图

    const int rough = ShaderManager::CheckSetShaderLocation(cloth_shader, "roughness");
    const int specStr = ShaderManager::CheckSetShaderLocation(cloth_shader, "specStrength");
    const int wrap = ShaderManager::CheckSetShaderLocation(cloth_shader, "wrapDiffuse");

    constexpr float clothRoughness = 0.80f;  // 越大越哑、highlight 越宽
    constexpr float clothSpec      = 0.22f;  // 高光强度：太大像塑料，太小没质感
    constexpr float clothWrap      = 0.25f;  // 漫反射包裹：增大可让暗面不至于太死
    SetShaderValue(cloth_shader, rough,   &clothRoughness, SHADER_UNIFORM_FLOAT);
    SetShaderValue(cloth_shader, specStr, &clothSpec,      SHADER_UNIFORM_FLOAT);
    SetShaderValue(cloth_shader, wrap,    &clothWrap,      SHADER_UNIFORM_FLOAT);
}

void HangingCloth::Update() {


}
