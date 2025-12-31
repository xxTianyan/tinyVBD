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
    const std::span<edge> edges;
    const std::span<triangle> tris;
    const std::span<tetrahedron> tets;
    const ForceElementAdjacencyInfo* adj = nullptr;
};

class World {
    // ! material and mesh must be one to one corresponded, which means
    //  m_materials.size() == mesh.size()
public:
    explicit World(Vec3  gravity);
    ~World()= default;

    void Add(MeshPtr m);

    void Remove();

    void Clear();

    static SimView MakeSimView(mesh_on_cpu& m);

    void InitStep();

    void ApplyGravity();

    void ChangeGravity(const Vec3& new_g){gravity = new_g;};

    MaterialID AddStVKMaterial(const StVkMaterial& m);

    std::vector<MeshPtr>  meshes;

    static bool RayNormal;

private:

    Vec3 gravity;
    size_t m_total_vertices{};
    size_t m_num_meshes{};

    // mesh -> material mapping (same length as meshes)
    std::vector<MaterialID> m_materials;

    std::vector<StVkMaterial> stvk_pool;


};

#endif //WORLD_H
