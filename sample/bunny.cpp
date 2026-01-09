//
// Created by tianyan on 1/9/26.
//

#include "bunny.h"

#include "Application.h"
#include "Builder.h"
#include "Scene.h"
#include "ShaderManager.h"
#include "VBDSolver.h"

void FallingBunny::CreateWorld(AppContext &ctx) {
    MModel model;
    Builder builder(model);
    m_bunny_id_ = builder.add_bunny();
    scene_ = std::make_unique<Scene>(std::move(model));
    solver_ = std::make_unique<VBDSolver>(&scene_->model_, 5, soft_bunny());
}

void FallingBunny::Step(const float dt) {
    solver_->Step(scene_->state_in(), scene_->state_out(), dt);
    scene_->SwapStates();
}


void FallingBunny::BindShaders(AppContext &ctx) {
    ctx.shader_manager->LoadShaderProgram("bunny", "../resources/shaders/bunny.vs", "../resources/shaders/bunny.fs");
    const auto bunny_shader = ctx.shader_manager->Get("bunny")->shader;
    ShaderManager::BindMatrices(bunny_shader);
    ShaderManager::SetCommonShaderParams(bunny_shader);
    bunny_shader.locs[SHADER_LOC_MAP_DIFFUSE] = GetShaderLocation(bunny_shader, "texture0");
    auto m_bunny_model = renderHelper_.GetRLModel(m_bunny_id_);
    m_bunny_model.materials[0].shader = bunny_shader;
    m_bunny_model.materials[0].maps[MATERIAL_MAP_ALBEDO].color = Color{230, 200, 160, 255};
}
