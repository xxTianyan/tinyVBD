//
// Created by 徐天焱 on 2025/11/11.
//

#include "Scene.h"

bool Scene::RayNormal = true;


/*MaterialID Scene::AddMaterial(const MMaterial& material) {
    const auto id = static_cast<MaterialID>(model_.materials.size());
    model_.materials.push_back(material);
    return id;
}

void Scene::BindMeshMaterial(const MeshID mesh, const MaterialID mat) {
    const auto mi = checked_index(mesh, model_.meshes, "mesh");
    const auto mid = checked_index(mat, model_.materials, "material");
    (void)mid;
    model_.mesh_to_material[mi] = mat;
}*/


/*void Scene::ApplyFixConsition(const MeshID mesh, const std::function<bool(const Vec3&)>& pred) {
    const auto mi = checked_index(mesh, model_.meshes, "mesh");

    auto& mm = model_.meshes[mi];
    auto& ms = state_.meshes[mi];

    // convention：fixed belongs to Model，pos/vel belong to State
    const size_t n = ms.pos.size();
    for (size_t i = 0; i < n; ++i) {
        if (pred(ms.pos[i])) {
            mm.fixed[i] = 1;
            mm.inv_mass[i] = 0.0f;

            // set dynamic value to 0
            ms.vel[i] = Vec3::Zero();
            ms.force[i] = Vec3::Zero();
            ms.prev_pos[i] = ms.pos[i];
            ms.inertia_pos[i] = ms.pos[i];
        }
    }
}*/

/*void Scene::InitStep() {
    ;
}*/

/*SimView Scene::MakeSimView(MeshID id) {
    auto& mm = model_.meshes.at(id);
    auto& ms = state_.meshes.at(id);

    auto mid = model_.mesh_to_material.at(id);
    return SimView{
        .pos = ms.pos,
        .prev_pos = ms.prev_pos,
        .inertia_pos = ms.inertia_pos,
        .vel = ms.vel,
        .accel = ms.force,
        .inv_mass = mm.inv_mass,
        .fixed = mm.fixed,
        .edges = mm.edges,
        .tris = mm.tris,
        .tets = mm.tets,
        .gravity = model_.gravity,
        .material_params = model_.materials.at(mid),
    };
}*/

























