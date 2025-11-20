//
// Created by 徐天焱 on 2025/11/11.
//

#ifndef WORLD_H
#define WORLD_H

#include "Mesh.h"

class World {

public:
    explicit World(Vec3  gravity);
    ~World()= default;

    void Step();

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
