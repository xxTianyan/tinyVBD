//
// Created by 徐天焱 on 2025/11/11.
//

#ifndef WORLD_H
#define WORLD_H

#include <span>
#include "MaterialParams.hpp"
#include "Mesh.h"

struct SimView {
    std::span<Vec3> pos;
    std::span<Vec3> prev_pos;
    std::span<Vec3> inertia_pos;
    std::span<Vec3> vel;
    std::span<Vec3> accel;
    std::span<Vec3> normal;
    std::span<float> inv_mass;
    const std::span<uint8_t> fixed;
    const std::span<edge> edges;
    const std::span<triangle> tris;
    const std::span<tetrahedron> tets;
    const Vec3& gravity;
    const MMaterial& material_params;
    const ForceElementAdjacencyInfo& adj;
};

class Scene {
public:
    explicit Scene(Vec3  gravity);
    ~Scene()= default;

    MeshID Add(MeshPtr m);

    void Remove();

    void Clear();

    void InitStep();

    void ChangeGravity(const Vec3& new_g){gravity = new_g;};

    void BindMeshMaterial(MeshID mesh, MaterialID mat);

    void ApplyFixConsition(MeshID mesh, const std::function<bool(const Vec3&)> &predicate);

    SimView MakeSimView(size_t mesh_id);

    MaterialID AddMaterial(MMaterial _params);

    // MeshID is mesh index
    std::vector<MeshPtr>  meshes;

    static bool RayNormal;

private:

    Vec3 gravity;
    size_t m_total_vertices{};
    size_t m_num_meshes{};

     // One material can bound to different meshes, but when this material is changed, all meshes bound with
     // this material will be influenced.
    std::vector<MMaterial> m_materials;
    // record mesh's material id, for meshes[mesh_id], its material is m_material[m_mesh_to_material[mesh_id]]
    std::vector<MaterialID> m_mesh_to_material;

};

#endif //WORLD_H
