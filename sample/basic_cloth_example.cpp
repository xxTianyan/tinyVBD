//
// Created by xumiz on 2026/1/1.
//

#include "basic_cloth_example.h"
#include "MeshBuilder.h"
#include "RenderHelper.hpp"

void HangingCloth::CreateWorld() {
    m_world->ChangeGravity(Vec3{0.0f, -9.81f, 0.0f});
    auto m = std::make_unique<mesh_on_cpu>();
    MeshBuilder::BuildCloth(m.get(), 2.0f, 2.0f, 10, 20, Vec3{0.0f, 3.0f, 0.0f}, ClothOrientation::Horizontal);
    if (!m_world) throw std::runtime_error("m_world is empty pointer");
    const auto mesh_id = m_world->Add(std::move(m));
    m_models = upload_all_models(*m_world);

    // create and bind material
    const auto material_id = m_world->AddMaterial(default_cloth());
    m_world->BindMeshMaterial(mesh_id, material_id);

    // fix boundary
    auto fix_left_z = [](const Vec3& pos) {
        if (std::abs(pos.z() - 1.0f) < 1e-5 ) return true;
        return false;
    };

    m_world->ApplyFixConsition(mesh_id, fix_left_z);
}

void HangingCloth::BindShaders() const {
    m_shader_manager->LoadShaderProgram("cloth", "../resources/shaders/cloth.vs", "../resources/shaders/cloth.fs");
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