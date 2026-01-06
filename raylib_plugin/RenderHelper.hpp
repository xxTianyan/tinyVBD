//
// Created by 徐天焱 on 2025/11/5.
//

#ifndef TAIYI_RENDERHELPER_H
#define TAIYI_RENDERHELPER_H


#include <algorithm>
#include <cstring>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <raylib.h>

#include "Model.h"
#include "Scene.h"

static void DrawAxisGizmo(const float length = 0.5f) {
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

static void FillVerticesXYZ(const State& state, const range particles, float* dst_xyz) {
    if (particles.end() > state.particle_pos.size()) {
        throw std::runtime_error("FillVerticesXYZ: particle range out of bounds");
    }

    for (size_t i = 0; i < particles.count; ++i) {
        const Vec3& p = state.particle_pos[particles.begin + i];
        dst_xyz[i * 3 + 0] = p.x();
        dst_xyz[i * 3 + 1] = p.y();
        dst_xyz[i * 3 + 2] = p.z();
    }
}

// compute real-time vertex normal and upload to gpu later
static void FillNormalsXYZ(const MModel& model,
                           const State& state,
                           const range& tri_range,
                           const range& particle_range,
                           float* dst_nxyz) {

    if (tri_range.end() > model.tris.size()) {
        throw std::runtime_error("FillNormalsXYZ: triangle range out of bounds");
    }
    if (particle_range.end() > state.particle_pos.size()) {
        throw std::runtime_error("FillNormalsXYZ: particle range out of bounds");
    }

    std::memset(dst_nxyz, 0, particle_range.count * 3 * sizeof(float));

    for (size_t t = 0; t < tri_range.count; ++t) {
        const triangle& tri = model.tris[tri_range.begin + t];

        const auto i0 = static_cast<size_t>(tri.vertices[0]);
        const auto i1 = static_cast<size_t>(tri.vertices[1]);
        const auto i2 = static_cast<size_t>(tri.vertices[2]);

        if (i0 < particle_range.begin || i1 < particle_range.begin || i2 < particle_range.begin) {
            throw std::runtime_error("FillNormalsXYZ: triangle references vertex before range start");
        }

        const size_t l0 = i0 - particle_range.begin;
        const size_t l1 = i1 - particle_range.begin;
        const size_t l2 = i2 - particle_range.begin;

        if (l0 >= particle_range.count || l1 >= particle_range.count || l2 >= particle_range.count) {
            throw std::runtime_error("FillNormalsXYZ: triangle vertex out of range");
        }

        const Vec3& p0 = state.particle_pos[i0];
        const Vec3& p1 = state.particle_pos[i1];
        const Vec3& p2 = state.particle_pos[i2];

        const Vec3 e1 = p1 - p0;
        const Vec3 e2 = p2 - p0;
        const Vec3 fn = e1.cross(e2); // 面法线（面积加权）

        dst_nxyz[l0 * 3 + 0] += fn.x();
        dst_nxyz[l0 * 3 + 1] += fn.y();
        dst_nxyz[l0 * 3 + 2] += fn.z();

        dst_nxyz[l1 * 3 + 0] += fn.x();
        dst_nxyz[l1 * 3 + 1] += fn.y();
        dst_nxyz[l1 * 3 + 2] += fn.z();

        dst_nxyz[l2 * 3 + 0] += fn.x();
        dst_nxyz[l2 * 3 + 1] += fn.y();
        dst_nxyz[l2 * 3 + 2] += fn.z();
    }

    // normalize per-vertex
    for (size_t i = 0; i < particle_range.count; ++i) {
        const float nx = dst_nxyz[i * 3 + 0];
        const float ny = dst_nxyz[i * 3 + 1];
        const float nz = dst_nxyz[i * 3 + 2];

        if (const float len2 = nx * nx + ny * ny + nz * nz; len2 > 1e-30f) {
            const float inv = 1.0f / std::sqrt(len2);
            dst_nxyz[i * 3 + 0] = nx * inv;
            dst_nxyz[i * 3 + 1] = ny * inv;
            dst_nxyz[i * 3 + 2] = nz * inv;
        } else {
            dst_nxyz[i * 3 + 0] = 0.0f;
            dst_nxyz[i * 3 + 1] = 1.0f;
            dst_nxyz[i * 3 + 2] = 0.0f;
        }
    }
}

// upload all mesh info into GPU Mesh + Model (dynamic)
static Model upload_model_from_cpu_mesh(const MModel& model, const State& state, const MeshInfo& info) {

    Mesh gmsh{};
    const size_t num_verts = info.particle.count;
    const size_t tri_count = info.tri.count;
    const size_t num_idx   = tri_count * 3;

    gmsh.vertexCount   = to_int_checked(num_verts, "vertexCount");
    gmsh.triangleCount = to_int_checked(tri_count, "triangleCount");

    // ---- vertices ----
    const size_t v_float_count = num_verts * 3;
    gmsh.vertices = static_cast<float*>(MemAlloc(v_float_count * sizeof(float)));
    if (!gmsh.vertices) throw std::runtime_error("MemAlloc failed: vertices");
    FillVerticesXYZ(state, info.particle, gmsh.vertices);

    // ---- indices (raylib Mesh.indices is unsigned short*) ----
    if (info.particle.count > static_cast<size_t>(std::numeric_limits<unsigned short>::max()) + 1ull) {
        throw std::runtime_error("Too many vertices for 16-bit indices (raylib Mesh.indices is unsigned short).");
    }

    static_assert(sizeof(unsigned short) == sizeof(uint16_t),
              "raylib Mesh.indices is expected to be 16-bit (unsigned short).");

    std::vector<uint16_t> idx;
    idx.reserve(num_idx);
    for (size_t t = 0; t < tri_count; ++t) {
        const triangle& tri = model.tris[info.tri.begin + t];
        for (int k = 0; k < 3; ++k) {
            const auto global = static_cast<size_t>(tri.vertices[k]);
            if (global < info.particle.begin || global >= info.particle.end()) {
                throw std::runtime_error("upload_model_from_cpu_mesh: triangle vertex out of range");
            }
            idx.push_back(static_cast<uint16_t>(global - info.particle.begin));
        }
    }

    gmsh.indices = static_cast<unsigned short*>(MemAlloc(num_idx * sizeof(unsigned short)));
    if (!gmsh.indices) throw std::runtime_error("MemAlloc failed: indices");
    std::memcpy(gmsh.indices, idx.data(), idx.size() * sizeof(uint16_t));

    // ---- normals (optional) ----
    if (Scene::RayNormal) {
        gmsh.normals = static_cast<float*>(MemAlloc(v_float_count * sizeof(float)));
        if (!gmsh.normals) throw std::runtime_error("MemAlloc failed: normals");
        FillNormalsXYZ(model, state, info.tri, info.particle, gmsh.normals);
    }

    UploadMesh(&gmsh, /*dynamic=*/true);

    // Attention：LoadModelFromMesh will take gmsh memory（raylib convention），
    // do not release gmsh.vertices/indices/normals maunally
    const Model rl_model = LoadModelFromMesh(gmsh);
    return rl_model;
}

// put every mesh to a model vector
static std::vector<Model> upload_all_models(const Scene& scene) {
    const auto& model = scene.model_;
    const State& state = scene.state_out_.particle_pos.empty() ? scene.state_in_ : scene.state_out_;

    std::vector<Model> out;
    out.reserve(model.mesh_infos.size());
    for (const auto& info : model.mesh_infos) {
        out.push_back(upload_model_from_cpu_mesh(model, state, info));
    }
    return out;
}

// update gpu mesh according to cpu mesh
static void UpdateModel(const std::vector<Model>& gpu_models, const Scene& scene) {
    const auto& cpu_model = scene.model_;
    const State& state = scene.state_out_.particle_pos.empty() ? scene.state_in_ : scene.state_out_;

    if (gpu_models.size() != cpu_model.mesh_infos.size()) {
        throw std::runtime_error("Model size mismatch");
    }

    std::vector<float> vertices; // 复用缓冲，避免每帧 malloc/free
    std::vector<float> normals;

    for (size_t i = 0; i < gpu_models.size(); ++i) {
        const MeshInfo& info = cpu_model.mesh_infos[i];

        const size_t float_count = info.particle.count * 3;

        // --- vertices ---
        vertices.resize(float_count);
        FillVerticesXYZ(state, info.particle, vertices.data());

        const size_t v_bytes = vertices.size() * sizeof(float);
        if (v_bytes > static_cast<size_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("Vertex buffer too large for int parameter");
        }

        // raylib: meshes[0] 是主 mesh
        const Mesh& gpu_mesh = gpu_models[i].meshes[0];

        // slot 0 = vertices
        UpdateMeshBuffer(gpu_mesh, 0, vertices.data(), static_cast<int>(v_bytes), 0);

        // --- normals (optional) ---
        if (Scene::RayNormal) {
            normals.resize(float_count);
            FillNormalsXYZ(cpu_model, state, info.tri, info.particle, normals.data());

            const size_t n_bytes = normals.size() * sizeof(float);
            if (n_bytes > static_cast<size_t>(std::numeric_limits<int>::max())) {
                throw std::runtime_error("Normal buffer too large for int parameter");
            }

            // slot 2 = normals
            UpdateMeshBuffer(gpu_mesh, 2, normals.data(), static_cast<int>(n_bytes), 0);
        }
    }
}



#endif //TAIYI_RENDERHELPER_H
