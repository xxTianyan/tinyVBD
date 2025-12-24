//
// Created by xumiz on 2025/12/24.
//

#include "Sample.h"
#include "MeshBuilder.h"
#include "RenderHelper.hpp"

/*
 * TODO: Finish Sample class and use new shader from gpt.
 */

HangingCloth::HangingCloth() {
    m_world = std::make_unique<World>(Vec3{0.0f, -9.81f, 0.0f});
    m_shader_manager = std::make_unique<ShaderManager>();
    m_solver = std::make_unique<VBDSolver>(10);
}

void HangingCloth::CreateWorld() {
    auto m = std::make_unique<mesh_on_cpu>();
    MeshBuilder::BuildCloth(m.get(), 1.0f, 2.0f, 10, 20, Vec3{0.0f, 2.0f, 0.0f});
    m_world->Add(std::move(m));
    m_models = upload_all_models(m_world.get());
}

void HangingCloth::BindShaders() {
    m_shader_manager->LoadShaderProgram("cloth", "../shaders/cloth.vs", "../shaders/cloth.fs");
    m_shader_manager->LoadShaderProgram("floor", "../shaders/floor.vs", "../shaders/floor.fs");
}

void HangingCloth::Update() {


}
