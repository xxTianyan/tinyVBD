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
    std::span<Vec3> force;
    std::span<Vec3> normal;
    std::span<Vec3> inv_m;
    /*const std::span<tetrahedron> tets;
    const std::span<triangle> tris;
    const std::span<edge> edges;*/
    const ForceElementAdjacencyInfo* adj = nullptr;
};

class World {

public:
    explicit World(Vec3  gravity);
    ~World()= default;

    static SimView MakeSimView(mesh_on_cpu& m);

    void Add(MeshPtr m);

    void Remove();

    void Clear();

    std::vector<MeshPtr>  meshes;

    static bool RayNormal;

private:

    Vec3 gravity;
    size_t m_total_vertices{};
    size_t m_num_meshes{};


};

#endif //WORLD_H
