//
// Created by xumiz on 2026/1/1.
//

#include "basic_cloth_example.h"
#include "Application.h"
#include "MeshBuilder.h"
#include "rlgl.h"
#include "Scene.h"
#include "VBDDynamics.h"

void BasicCloth::CreateWorld([[maybe_unused]]AppContext &ctx) {
    scene_ = std::make_unique<Scene>(Vec3{0.0f,-9.81f, 0.0f});
    solver_ = std::make_unique<VBDSolver>(10);
    auto m = std::make_unique<mesh_on_cpu>();
    MeshBuilder::BuildCloth(m.get(), 2.0f, 2.0f, 10, 20, Vec3{0.0f, 3.0f, 0.0f}, ClothOrientation::Horizontal);
    const auto mesh_id = scene_->Add(std::move(m));

    // create and bind material
    const auto material_id = scene_->AddMaterial(default_cloth());
    scene_->BindMeshMaterial(mesh_id, material_id);

    auto fix_left_z = [](const Vec3& pos) {
        if (std::abs(pos.z() - 1.0f) < 1e-5 ) return true;
        return false;
    };

    scene_->ApplyFixConsition(mesh_id, fix_left_z);

}

void BasicCloth::Render(AppContext &ctx) {
    BeginMode3D(ctx.orbitCam->camera);

    // floor
    if (IsModelValid_(floor_)) {
        DrawModel(floor_, Vector3{0,0,0}, 1.0f, WHITE);
    }

    DrawAxisGizmo();

    // scene models
    for (auto& m : models_) {
        if (!IsModelValid_(m)) continue;
        rlDisableBackfaceCulling();
        DrawModel(m, Vector3{0,0,0}, 1.0f, WHITE);
        rlEnableBackfaceCulling();
    }

    EndMode3D();
}


void BasicCloth::BindShaders(AppContext &ctx) {
    ctx.shader_manager->LoadShaderProgram("cloth", "../resources/shaders/cloth.vs", "../resources/shaders/cloth.fs");
    const auto cloth_shader = ctx.shader_manager->Get("cloth")->shader;
    ShaderManager::BindMatrices(cloth_shader);
    ShaderManager::SetCommonShaderParams(cloth_shader);
    cloth_shader.locs[SHADER_LOC_MAP_DIFFUSE] = GetShaderLocation(cloth_shader, "texture0");
    auto m_cloth = models_[0];
    m_cloth.materials[0].shader = cloth_shader;
    m_cloth.materials[0].maps[MATERIAL_MAP_ALBEDO].color = Color{230, 200, 160, 255};
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

void BasicCloth::BuildRenderResources() {
    models_ = upload_all_models(*scene_);
}

void BasicCloth::Step(const float dt) {
    scene_->InitStep();
    auto& meshes = scene_->meshes;
    for (size_t mesh_id = 0; mesh_id < meshes.size(); mesh_id++) {
        SimView view = scene_->MakeSimView(mesh_id);
        // make inertia step
        VBDSolver::forward_step(view, dt);
        // iter newton step
        for (size_t iter = 0; iter < 5; iter++)
            VBDSolver::solve(view, dt);
        VBDSolver::update_velocity(view, dt);
    }
}

