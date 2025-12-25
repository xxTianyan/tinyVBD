//
// Created by 徐天焱 on 2025/11/5.
//

#ifndef RENDERHELPER_H
#define RENDERHELPER_H

#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include <raylib.h>

#include "Mesh.h"
#include "World.h"

static void DrawAxisGizmo(float length = 0.5f) {
    constexpr Vector3 origin = { 0.0f, 0.0f, 0.0f };
    const Vector3 x = { length, 0.0f, 0.0f };
    const Vector3 y = { 0.0f, length, 0.0f };
    const Vector3 z = { 0.0f, 0.0f, length };

    DrawLine3D(origin, x, RED);
    DrawLine3D(origin, y, GREEN);
    DrawLine3D(origin, z, BLUE);

    const float radius = length * 0.03f;
    DrawSphereWires(x, radius, 1, 6, RED);
    DrawSphereWires(y, radius, 1, 6, GREEN);
    DrawSphereWires(z, radius, 1, 6, BLUE);
}

// turn mesh_on_cpu into GPU Mesh + Model (dynamic)
static Model upload_model_from_cpu_mesh(mesh_on_cpu* M) {
    Mesh gmsh = {0};
    gmsh.vertexCount   = static_cast<int>(M->size());
    gmsh.triangleCount = static_cast<int>(M->m_surface_tris.size() / 3);

    const std::vector<float> vertices = assemble_vertices(M);
    // copy vertices
    gmsh.vertices = static_cast<float *>(MemAlloc(vertices.size() * sizeof(float)));
    std::memcpy(gmsh.vertices, vertices.data(), vertices.size() * sizeof(float));
    // copy indices
    gmsh.indices = static_cast<unsigned short *>(MemAlloc(M->m_surface_tris.size() * sizeof(unsigned short)));
    std::memcpy(gmsh.indices, M->m_surface_tris.data(), M->m_surface_tris.size() * sizeof(unsigned short));

    if (World::RayNormal) {
        const std::vector<float> normals  = ComputeNormal(M);
        gmsh.normals = static_cast<float *>(MemAlloc(normals.size() * sizeof(float)));
        std::memcpy(gmsh.normals, normals.data(), normals.size() * sizeof(float));
    }

    UploadMesh(&gmsh, /*dynamic=*/true);

    // raylib will hold the mesh
    const Model model = LoadModelFromMesh(gmsh);
    return model;
}

// put every mesh to a model vector
static std::vector<Model> upload_all_models(const World& world) {
    std::vector<Model> models;
    models.reserve(world.meshes.size());
    for (auto& up : world.meshes) {
        models.push_back(upload_model_from_cpu_mesh(up.get()));
    }
    return models;
}

// update gpu mesh according to cpu mesh
static void UpdateModel(const std::vector<Model> &models, const std::vector<MeshPtr>& Meshes) {
    if (models.size() != Meshes.size()) throw std::runtime_error("Model size mismatch");
    std::vector<float> vertices;
    std::vector<float> normals;
    for (size_t i = 0; i < models.size(); ++i) {
        const Mesh& old_mesh = models[i].meshes[0];
        vertices = assemble_vertices(Meshes[i].get());
        normals = ComputeNormal(Meshes[i].get());
        const size_t dataSize = vertices.size() * sizeof(float);
        if (dataSize > static_cast<size_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("Data too large for int parameter");
        }
        UpdateMeshBuffer(old_mesh, 0, vertices.data(), static_cast<int>(dataSize),0);
        UpdateMeshBuffer(old_mesh, 2, normals.data(), static_cast<int>(dataSize), 0);
    }
}

#endif //RENDERHELPER_H
