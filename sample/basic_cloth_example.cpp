//
// Created by xumiz on 2026/1/1.
//

#include "basic_cloth_example.h"
#include "Application.h"
#include "Builder.h"
#include "rlgl.h"
#include "Scene.h"
#include "VBDSolver.h"

void BasicCloth::CreateWorld([[maybe_unused]]AppContext &ctx) {
    MModel model;
    Builder builder(model);
    m_cloth_id_ = builder.add_cloth(2.0f, 3.0f, 16, 24, Vec3{0.0f, 4.0f, 0.0f});
    scene_ = std::make_unique<Scene>(std::move(model));
    solver_ = std::make_unique<VBDSolver>(&scene_->model_, 5);
}

void BasicCloth::Render(AppContext &ctx) {

    BeginMode3D(ctx.orbitCam->camera);

    // floor
    if (RenderHelper::IsModelValid(floor_)) {
        DrawModel(floor_, Vector3{0,0,0}, 1.0f, WHITE);
    }

    // scene models
    rlDisableBackfaceCulling();
    renderHelper_.Draw(ctx.is_wire_mode);
    rlEnableBackfaceCulling();

    EndMode3D();
}


void BasicCloth::BindShaders(AppContext &ctx) {
    ctx.shader_manager->LoadShaderProgram("cloth", "../resources/shaders/cloth.vs", "../resources/shaders/cloth.fs");
    const auto cloth_shader = ctx.shader_manager->Get("cloth")->shader;
    ShaderManager::BindMatrices(cloth_shader);
    ShaderManager::SetCommonShaderParams(cloth_shader);
    cloth_shader.locs[SHADER_LOC_MAP_DIFFUSE] = GetShaderLocation(cloth_shader, "texture0");
    auto m_cloth_model = renderHelper_.GetRLModel(m_cloth_id_);
    m_cloth_model.materials[0].shader = cloth_shader;
    m_cloth_model.materials[0].maps[MATERIAL_MAP_ALBEDO].color = Color{230, 200, 160, 255};
    // m_cloth.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = clothAlbedoTex;          // if texture

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


void BasicCloth::Step(const float dt) {
    scene_->InitStep();
    solver_->Step(scene_->state_in(), scene_->state_out(), dt);
    scene_->SwapStates();
}

