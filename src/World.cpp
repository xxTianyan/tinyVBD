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


void World::Step() {
    for (const auto& m : meshes) {
        for (size_t i = 0; i < m->size(); ++i) {
            m->px_pred[i] = m->px[i] + (gravity[0] + m->fx[i]) * 1.0f/60; // at random
        }
    }
}




