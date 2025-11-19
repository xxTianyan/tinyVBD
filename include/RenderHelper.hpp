//
// Created by 徐天焱 on 2025/11/5.
//

#ifndef HELPER_H
#define HELPER_H

#include <iostream>
#include <raylib.h>
#include "raymath.h"
#include "Mesh.h"
#include"World.h"

template<typename T>
void print_vector(const std::vector<T>& vec) {
    for (const auto& elem : vec) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}

struct OrbitCtrl {
    float yaw;        // 水平角 (弧度)
    float pitch;      // 俯仰角 (弧度)
    float radius;     // 目标点到相机的距离
    float rotSens;    // 旋转灵敏度
    float zoomSens;   // 缩放灵敏度
    float panSens;    // 平移灵敏度(像素->世界)
};

inline void ClampPitch(float *pitch) {
    constexpr float lim = DEG2RAD*89.0f;                 // 避免万向节锁
    if (*pitch >  lim) *pitch =  lim;
    if (*pitch < -lim) *pitch = -lim;
}

inline Vector3 SphericalToCartesian(const float r, const float yaw, const float pitch) {
    const float cp = cosf(pitch);
    const float sp = sinf(pitch);
    const float cy = cosf(yaw);
    const float sy = sinf(yaw);
    return (Vector3){ r*cp*cy, r*sp, r*cp*sy };
}

inline void ReframeToModel(Camera3D *cam,
                           OrbitCtrl *orbit,
                           const std::vector<Model> &models,
                           const float margin)
{
    if (!cam || !orbit || models.empty()) return;

    BoundingBox box0 = GetModelBoundingBox(models[0]);
    Vector3 min = box0.min;
    Vector3 max = box0.max;

    for (size_t i = 1; i < models.size(); ++i) {
        BoundingBox bi = GetModelBoundingBox(models[i]);
        min.x = fminf(min.x, bi.min.x);
        min.y = fminf(min.y, bi.min.y);
        min.z = fminf(min.z, bi.min.z);

        max.x = fmaxf(max.x, bi.max.x);
        max.y = fmaxf(max.y, bi.max.y);
        max.z = fmaxf(max.z, bi.max.z);
    }

    const Vector3 center = {
        (min.x + max.x) * 0.5f,
        (min.y + max.y) * 0.5f,
        (min.z + max.z) * 0.5f
    };

    const Vector3 diag = Vector3Subtract(max, min);
    const float radiusBox = 0.5f * Vector3Length(diag);
    const float fitDist   = radiusBox / tanf(DEG2RAD * cam->fovy * 0.5f);

    const float m = (margin > 1.0f ? margin : 1.15f);
    orbit->radius = fitDist * m;

    cam->target = center;

    const Vector3 offset = SphericalToCartesian(orbit->radius, orbit->yaw, orbit->pitch);
    cam->position = Vector3Add(cam->target, offset);
    cam->up       = (Vector3){ 0.0f, 1.0f, 0.0f };
}


// 把一个 mesh_on_cpu 变成 GPU Mesh + Model（动态可更新）
static Model upload_model_from_cpu_mesh(mesh_on_cpu* M)
{
    Mesh gmsh = {0};
    gmsh.vertexCount   = static_cast<int>(M->size());
    gmsh.triangleCount = static_cast<int>(M->m_surface_tris_local.size() / 3);

    const std::vector<float> vertices = assemble_vertices(M);
    // copy vertices
    gmsh.vertices = static_cast<float *>(MemAlloc(vertices.size() * sizeof(float)));
    std::memcpy(gmsh.vertices, vertices.data(), vertices.size() * sizeof(float));
    // copy indices
    gmsh.indices = static_cast<unsigned short *>(MemAlloc(M->m_surface_tris_local.size() * sizeof(unsigned short)));
    std::memcpy(gmsh.indices, M->m_surface_tris_local.data(), M->m_surface_tris_local.size() * sizeof(unsigned short));

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


// 把 World 里的每个 mesh_on_cpu 都上传成一个 Model
static std::vector<Model> upload_all_models(const World& world)
{
    std::vector<Model> models;
    models.reserve(world.meshes.size());
    for (auto& up : world.meshes) {
        models.push_back(upload_model_from_cpu_mesh(up.get()));
    }
    return models;
}


#endif //HELPER_H
