//
// Created by tianyan on 1/9/26.
//

#include "bunny.h"

#include <iostream>

#include "Application.h"
#include "Builder.h"
#include "Scene.h"
#include "ShaderManager.h"
#include "VBDSolver.h"

void print_quality_histogram(const std::vector<tetrahedron>& tets) {
    int buckets[10] = {0}; // 0.0-0.1, 0.1-0.2, ... 0.9-1.0
    int inverted_count = 0;

    for (const auto& t : tets) {
        float q = t.qualityStats.meanRatio;
        if (t.qualityStats.scaledJacobian < 0) {
            inverted_count++;
        }

        if (q >= 0.0f && q <= 1.0f) {
            int idx = static_cast<int>(q * 10);
            if (idx == 10) idx = 9; // Handle exact 1.0
            buckets[idx]++;
        }
    }

    std::cout << "\n=== Mesh Quality Report (Mean Ratio) ===" << std::endl;
    std::cout << "Inverted Elements (Jacobian < 0): " << inverted_count << " !!!" << std::endl;

    int max_count = 0;
    for(int c : buckets) max_count = std::max(max_count, c);

    for (int i = 0; i < 10; ++i) {
        printf("[%0.1f - %0.1f]: %6d | ", i/10.0f, (i+1)/10.0f, buckets[i]);

        // 简单的 ASCII 条形图
        int bar_len = (int)((float)buckets[i] / max_count * 20.0f);
        for(int k=0; k<bar_len; ++k) std::cout << "#";
        std::cout << std::endl;
    }
    std::cout << "========================================" << std::endl;
}

void FallingBunny::CreateWorld(AppContext &ctx) {
    MModel model;
    Builder builder(model);
    m_bunny_id_ = builder.add_sphere(1.0f, 10, Vec3{0.0f,30.0f,0.0f}, 3.f, "sphere");
    // m_bunny_id_ = builder.add_bunny();
    scene_ = std::make_unique<Scene>(std::move(model));
    solver_ = std::make_unique<VBDSolver>(&scene_->model_, 10, soft_bunny());
    print_quality_histogram(scene_->model_.tets);
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
