//
// Created by 徐天焱 on 2025/11/11.
//

#ifndef WORLD_H
#define WORLD_H

#include <span>
#include "Mesh.h"

struct SimView {
    std::span<Vec3> pos;
    std::span<Vec3> pred_pos;
    std::span<Vec3> inertia_pos;
    std::span<Vec3> vel;
    std::span<Vec3> accel;
    std::span<Vec3> normal;
    std::span<Vec3> inv_m;
    const std::span<edge> edges;
    const std::span<triangle> tris;
    const std::span<tetrahedron> tets;
    const ForceElementAdjacencyInfo* adj = nullptr;
};

class World {

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

    std::vector<MeshPtr>  meshes;

    static bool RayNormal;

private:

    Vec3 gravity;
    size_t m_total_vertices{};
    size_t m_num_meshes{};


};

#endif //WORLD_H
