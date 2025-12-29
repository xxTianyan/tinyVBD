//
// Created by 徐天焱 on 2025/11/11.
//

#include "World.h"

bool World::RayNormal = true;

World::World(Vec3  gravity) : gravity(std::move(gravity)) {}

void World::Add(MeshPtr m) {
    m->base_offset = m_total_vertices;
    m_total_vertices += m->size();
    BuildAdjacency(*m);
    meshes.push_back(std::move(m));
}

void World::Clear() {
    ;
}

void World::Remove() {
    ;
}

SimView World::MakeSimView(mesh_on_cpu& m) {
    return SimView{
        .pos = m.p,
        .pred_pos = m.p_pred,
        .inertia_pos = m.p_inertia,
        .vel = m.v,
        .accel = m.accel,
        .normal = m.n,
        .inv_m = m.inv_m,
        .adj = &m.adjacencyInfo
    };
}

void World::ApplyGravity() {
    for (auto& m : meshes) {
        auto& accel = m->accel;
        std::ranges::fill(accel.begin(), accel.end(), gravity);
    }
}

void World::InitStep() {
    ApplyGravity();
}

























