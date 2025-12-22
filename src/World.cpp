//
// Created by 徐天焱 on 2025/11/11.
//

#include "World.h"

bool World::RayNormal = true;

World::World(Vec3  gravity) : gravity(std::move(gravity)) {}


void World::Add(MeshPtr m) {
    m->base_offest = m_total_vertices;
    m_total_vertices += m->size();
    meshes.push_back(std::move(m));
}

void World::Clear() {
    meshes.clear();
}

void World::Remove() {


}


void World::Step(const float dt) {
    for (const auto& m : meshes) {
        for (size_t i = 0; i < m->size(); ++i) {
            m->py[i] +=  0.1f * dt; // at random
        }
    }
}




