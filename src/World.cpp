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
        .pos = m.pos,
        .pred_pos = m.pred_pos,
        .inertia_pos = m.inertia_pos,
        .vel = m.v,
        .accel = m.accel,
        .normal = m.n,
        .inv_m = m.inv_m,
        .edges = m.m_edges,
        .tris = m.m_tris,
        .tets = m.m_tets,
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

























