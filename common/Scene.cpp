//
// Created by 徐天焱 on 2025/11/11.
//

#include "Scene.h"

bool Scene::RayNormal = true;

Scene::Scene(Vec3  gravity) : gravity(std::move(gravity)) {}

MeshID Scene::Add(MeshPtr m) {
    m->base_offset = m_total_vertices;
    m_total_vertices += m->size();
    InitMesh(*m);
    meshes.push_back(std::move(m));
    return static_cast<MeshID>(meshes.size()-1);
}

void Scene::Clear() {
    ;
}

void Scene::Remove() {
    ;
}

void Scene::ApplyFixConsition(MeshID _mid, const std::function<bool(const Vec3&)> &predicate) {
    const auto& m = meshes[_mid];
    for (size_t i = 0; i < m->size(); ++i) {
        if (const auto& pos = m->pos[i]; predicate(pos))
            m->fixed[i] = 1;
    }
}

SimView Scene::MakeSimView(const size_t mesh_id) {
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
        .fixed = m.fixed,
        .edges = m.m_edges,
        .tris = m.m_tris,
        .tets = m.m_tets,
        .gravity = gravity,  // gravity is owned by world
        .material_params = m_materials[m_mesh_to_material[mesh_id]],
        .adj = m.adjacencyInfo
    };
}

void Scene::InitStep() {
    for (auto& m : meshes)
        std::fill(m->accel.begin(), m->accel.end(),Vec3::Zero());
}

MaterialID Scene::AddMaterial(MMaterial _params) {
    // valid prams check
    if (!(_params.E() > 0.0f)) throw std::runtime_error("Young's module must be > 0");
    if (!(_params.nu() > -0.99f && _params.nu() < 0.49f)) throw std::runtime_error("Poisson's ratio out of range");
    // push to my materials
    m_materials.emplace_back(_params);
    return static_cast<uint32_t>(m_materials.size()) - 1;
}

void Scene::BindMeshMaterial(const MeshID mesh, const MaterialID mat) {
    if (mesh >= meshes.size()) throw std::out_of_range("Mesh ID is out of range");
    if (mat >= m_materials.size()) throw std::out_of_range("Material ID is out of range");
    m_mesh_to_material.push_back(mat);
}

























