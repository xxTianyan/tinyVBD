//
// Created by tianyan on 1/6/26.
//

#include "RenderHelper.h"


void RenderHelper::BindModel(const MModel& model) {
    model_ = &model;
    // 这里不重建，避免调用端 Bind 后还没准备 state 就触发构建
    // 真正构建在 Update() 中按版本判断触发
}

void RenderHelper::CheckBound() const {
    if (!model_) throw std::runtime_error("RenderHelper: model not bound. Call BindModel(model) first.");
}

void RenderHelper::Shutdown() {
    for (auto& rm : meshes_) {
        if (rm.valid) {
            UnloadModel(rm.model);
            rm.valid = false;
        }
    }
    meshes_.clear();
    ready_ = false;
    built_topology_version_ = 0;
}

void RenderHelper::Update(const State& state) {
    CheckBound();

    // 1) 基础一致性检查：state 必须覆盖 model.num_particles
    if (model_->num_particles != state.particle_pos.size()) {
        throw std::runtime_error("RenderHelper::Update: model.num_particles != state.particle_pos.size()");
    }

    // 2) 拓扑版本变化：重建
    if (!ready_ || built_topology_version_ != model_->topology_version) {
        Rebuild();
        built_topology_version_ = model_->topology_version;
        ready_ = true;
    }

    // 3) 更新 VBO（positions + normals）
    UpdateDynamic(state);
}

void RenderHelper::Rebuild() {
    CheckBound();

    // 先释放旧资源（安全）
    Shutdown();

    meshes_.reserve(model_->mesh_infos.size());

    for (const MeshInfo& info_src : model_->mesh_infos) {
        RenderMesh rm{};
        rm.info = info_src; // 值拷贝 ranges（避免 mesh_infos realloc 引用失效风险）

        const size_t particle_begin = rm.info.particle.begin;
        const size_t particle_count = rm.info.particle.count;
        const size_t tri_count      = rm.info.tri.count;

        if (particle_count == 0 || tri_count == 0) {
            rm.valid = false;
            meshes_.push_back(rm);
            continue;
        }

        // raylib Mesh indices 是 u16：单 mesh 顶点数限制
        if (particle_count > static_cast<size_t>(std::numeric_limits<unsigned short>::max()) + 1ull) {
            throw std::runtime_error("RenderHelper::Rebuild: particle_count > 65536, raylib u16 indices not supported. Split mesh.");
        }

        // range 合法性
        if (particle_begin + particle_count > model_->num_particles) {
            throw std::runtime_error("RenderHelper::Rebuild: particle range out of model.num_particles");
        }
        if (rm.info.tri.end() > model_->tris.size()) {
            throw std::runtime_error("RenderHelper::Rebuild: tri range out of model.tris");
        }

        // -------- build CPU mesh arrays (MemAlloc, owned by raylib) --------
        Mesh mesh{};
        mesh.vertexCount   = static_cast<int>(particle_count);
        mesh.triangleCount = static_cast<int>(tri_count);

        mesh.vertices = (float*)MemAlloc(particle_count * 3 * sizeof(float));
        mesh.normals  = (float*)MemAlloc(particle_count * 3 * sizeof(float));
        mesh.indices  = (unsigned short*)MemAlloc(tri_count * 3 * sizeof(unsigned short));

        if (!mesh.vertices || !mesh.normals || !mesh.indices) {
            throw std::runtime_error("RenderHelper::Rebuild: MemAlloc failed");
        }

        // 这里用 particle_pos0 初始化也可以；但更稳妥做法是：首次重建时也由 Update(state) 负责写
        // 为避免多做一次拷贝，我们先填充一个合理初值（零），随后 UpdateDynamic 会覆盖
        std::memset(mesh.vertices, 0, particle_count * 3 * sizeof(float));
        std::memset(mesh.normals,  0, particle_count * 3 * sizeof(float));

        BuildIndicesU16(*model_, rm.info.tri, particle_begin, particle_count, (unsigned short*)mesh.indices);

        // Upload as dynamic: positions/normals 每帧更新
        UploadMesh(&mesh, true);

        rm.model = LoadModelFromMesh(mesh);
        rm.model.materials[0].maps[MATERIAL_MAP_DIFFUSE].color = WHITE;

        rm.valid = true;
        meshes_.push_back(rm);
    }

    ready_ = true;
}

void RenderHelper::UpdateDynamic(const State& state) {
    CheckBound();
    if (!ready_) return;

    for (auto& rm : meshes_) {
        if (!rm.valid) continue;

        const size_t particle_begin = rm.info.particle.begin;
        const size_t particle_count = rm.info.particle.count;

        Mesh& mesh = rm.model.meshes[0];

        if (static_cast<size_t>(mesh.vertexCount) != particle_count) {
            throw std::runtime_error("RenderHelper::UpdateDynamic: mesh.vertexCount mismatch (topology out of sync)");
        }
        if (!mesh.vertices || !mesh.normals) {
            throw std::runtime_error("RenderHelper::UpdateDynamic: mesh CPU arrays are null");
        }

        FillPositionsXYZ(state, particle_begin, particle_count, mesh.vertices);
        ComputeNormalsXYZ(*model_, state, rm.info.tri, particle_begin, particle_count, mesh.normals);

        // raylib buffer slots: 0=positions, 2=normals
        UpdateMeshBuffer(mesh, 0, mesh.vertices, mesh.vertexCount * 3 * (int)sizeof(float), 0);
        UpdateMeshBuffer(mesh, 2, mesh.normals,  mesh.vertexCount * 3 * (int)sizeof(float), 0);
    }
}

void RenderHelper::Draw() const {
    if (!ready_) return;

    for (const auto& rm : meshes_) {
        if (!rm.valid) continue;
        DrawModel(rm.model, Vector3{0.0f, 0.0f, 0.0f}, 1.0f, WHITE);
    }
}

// ------------------------ math helpers ------------------------

void RenderHelper::FillPositionsXYZ(const State& state,
                                    const size_t particle_begin,
                                    const size_t particle_count,
                                    float* dst_xyz) {
    if (!dst_xyz) throw std::runtime_error("FillPositionsXYZ: dst_xyz null");
    if (particle_begin + particle_count > state.particle_pos.size()) {
        throw std::runtime_error("FillPositionsXYZ: particle range out of state.particle_pos");
    }

    for (size_t i = 0; i < particle_count; ++i) {
        const Vec3& p = state.particle_pos[particle_begin + i];
        dst_xyz[i * 3 + 0] = p.x();
        dst_xyz[i * 3 + 1] = p.y();
        dst_xyz[i * 3 + 2] = p.z();
    }
}

void RenderHelper::BuildIndicesU16(const MModel& model,
                                  const range tri_range,
                                  const size_t particle_begin,
                                  const size_t particle_count,
                                  unsigned short* dst_indices) {
    if (!dst_indices) throw std::runtime_error("BuildIndicesU16: dst null");
    if (tri_range.end() > model.tris.size()) {
        throw std::runtime_error("BuildIndicesU16: tri range out of model.tris");
    }

    for (size_t t = 0; t < tri_range.count; ++t) {
        const triangle& tri = model.tris[tri_range.begin + t];

        const VertexID g0 = tri.vertices[0];
        const VertexID g1 = tri.vertices[1];
        const VertexID g2 = tri.vertices[2];

        if (g0 < particle_begin || g1 < particle_begin || g2 < particle_begin) {
            throw std::runtime_error("BuildIndicesU16: triangle references vertex before particle_begin");
        }

        const size_t i0 = static_cast<size_t>(g0 - particle_begin);
        const size_t i1 = static_cast<size_t>(g1 - particle_begin);
        const size_t i2 = static_cast<size_t>(g2 - particle_begin);

        if (i0 >= particle_count || i1 >= particle_count || i2 >= particle_count) {
            throw std::runtime_error("BuildIndicesU16: triangle references vertex out of this mesh particle range");
        }

        dst_indices[t * 3 + 0] = static_cast<unsigned short>(i0);
        dst_indices[t * 3 + 1] = static_cast<unsigned short>(i1);
        dst_indices[t * 3 + 2] = static_cast<unsigned short>(i2);
    }
}

void RenderHelper::ComputeNormalsXYZ(const MModel& model,
                                     const State& state,
                                     const range tri_range,
                                     const size_t particle_begin,
                                     const size_t particle_count,
                                     float* dst_nxyz) {
    if (!dst_nxyz) throw std::runtime_error("ComputeNormalsXYZ: dst null");
    if (particle_begin + particle_count > state.particle_pos.size()) {
        throw std::runtime_error("ComputeNormalsXYZ: particle range out of state.particle_pos");
    }
    if (tri_range.end() > model.tris.size()) {
        throw std::runtime_error("ComputeNormalsXYZ: tri range out of model.tris");
    }

    std::memset(dst_nxyz, 0, particle_count * 3 * sizeof(float));

    constexpr float eps = 1e-12f;

    for (size_t t = 0; t < tri_range.count; ++t) {
        const triangle& tri = model.tris[tri_range.begin + t];

        const VertexID g0 = tri.vertices[0];
        const VertexID g1 = tri.vertices[1];
        const VertexID g2 = tri.vertices[2];

        // 映射到 mesh 局部 index
        if (g0 < particle_begin || g1 < particle_begin || g2 < particle_begin) continue;

        const size_t i0 = static_cast<size_t>(g0 - particle_begin);
        const size_t i1 = static_cast<size_t>(g1 - particle_begin);
        const size_t i2 = static_cast<size_t>(g2 - particle_begin);

        if (i0 >= particle_count || i1 >= particle_count || i2 >= particle_count) continue;

        const Vec3& p0 = state.particle_pos[particle_begin + i0];
        const Vec3& p1 = state.particle_pos[particle_begin + i1];
        const Vec3& p2 = state.particle_pos[particle_begin + i2];

        const Vec3 n = (p1 - p0).cross(p2 - p0); // area-weighted
        const float nx = n.x(), ny = n.y(), nz = n.z();

        dst_nxyz[i0 * 3 + 0] += nx; dst_nxyz[i0 * 3 + 1] += ny; dst_nxyz[i0 * 3 + 2] += nz;
        dst_nxyz[i1 * 3 + 0] += nx; dst_nxyz[i1 * 3 + 1] += ny; dst_nxyz[i1 * 3 + 2] += nz;
        dst_nxyz[i2 * 3 + 0] += nx; dst_nxyz[i2 * 3 + 1] += ny; dst_nxyz[i2 * 3 + 2] += nz;
    }

    // normalize
    for (size_t i = 0; i < particle_count; ++i) {
        const float x = dst_nxyz[i * 3 + 0];
        const float y = dst_nxyz[i * 3 + 1];
        const float z = dst_nxyz[i * 3 + 2];
        const float len2 = x * x + y * y + z * z;

        if (len2 > eps) {
            const float inv = 1.0f / std::sqrt(len2);
            dst_nxyz[i * 3 + 0] = x * inv;
            dst_nxyz[i * 3 + 1] = y * inv;
            dst_nxyz[i * 3 + 2] = z * inv;
        } else {
            // fallback
            dst_nxyz[i * 3 + 0] = 0.0f;
            dst_nxyz[i * 3 + 1] = 1.0f;
            dst_nxyz[i * 3 + 2] = 0.0f;
        }
    }
}

