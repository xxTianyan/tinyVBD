//
// Created by 徐天焱 on 2025/11/11.
//

#include "World.h"

bool World::RayNormal = true;

World::World(Vec3  gravity) : gravity(std::move(gravity)) {}

MeshID World::Add(MeshPtr m) {
    m->base_offset = m_total_vertices;
    m_total_vertices += m->size();
    InitMesh(*m);
    meshes.push_back(std::move(m));
    return static_cast<MeshID>(meshes.size()-1);
}

void World::Clear() {
    ;
}

void World::Remove() {
    ;
}

SimView World::MakeSimView(const size_t mesh_id) {
    if (mesh_id >= meshes.size()) throw std::out_of_range("Mesh ID is out of range");
    auto& m = *meshes[mesh_id];
    return SimView{
        .pos = m.pos,
        .prev_pos = m.prev_pos,
        .inertia_pos = m.inertia_pos,
        .vel = m.vel,
        .accel = m.accel,
        .normal = m.n,
        .inv_mass = m.inv_mass,
        .edges = m.m_edges,
        .tris = m.m_tris,
        .tets = m.m_tets,
        .gravity = gravity,  // gravity is owned by world
        .material_params = m_materials[m_mesh_to_material[mesh_id]],
        .adj = m.adjacencyInfo
    };
}

void World::ApplyGravity() {
    for (auto& m : meshes) {
        auto& accel = m->accel;
        std::ranges::fill(accel.begin(), accel.end(), gravity);
    }
}

void World::InitStep() {
    // ApplyGravity();
}

MaterialID World::AddMaterial(MMaterial _params) {
    // valid prams check
    if (!(_params.E() > 0.0f)) throw std::runtime_error("Young's module must be > 0");
    if (!(_params.nu() > -0.99f && _params.nu() < 0.49f)) throw std::runtime_error("Poisson's ratio out of range");
    // push to my materials
    m_materials.emplace_back(_params);
    return static_cast<uint32_t>(m_materials.size()) - 1;
}

void World::BindMeshMaterial(const MeshID mesh, const MaterialID mat) {
    if (mesh >= meshes.size()) throw std::out_of_range("Mesh ID is out of range");
    if (mat >= m_materials.size()) throw std::out_of_range("Material ID is out of range");
    m_mesh_to_material.push_back(mat);
}

























