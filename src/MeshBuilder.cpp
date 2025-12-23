//
// Created by tianyan on 12/23/25.
//

#include "MeshBuilder.h"
#include <cmath>

void MeshBuilder::PrepareMesh(mesh_on_cpu* mesh, const size_t num_nodes) {
    if (!mesh) return;
    mesh->clear_topology(); // 清空旧的索引数据
    mesh->resize(num_nodes); // 调整 SoA 数组大小
}

void MeshBuilder::BuildCloth(mesh_on_cpu* mesh, const float width, const float height, const int resX, const int resY, const Vec3& center) {
    const size_t num_nodes = (resX + 1) * (resY + 1);
    PrepareMesh(mesh, num_nodes);

    const float dx = width / resX;
    const float dy = height / resY;
    const Vec3 start_pos = center - Vec3(width * 0.5f, height * 0.5f, 0.0f);

    // 1. 生成顶点
    for (int j = 0; j <= resY; ++j) {
        for (int i = 0; i <= resX; ++i) {
            const int idx = j * (resX + 1) + i;
            mesh->p[idx] = start_pos + Vec3(i * dx, j * dy, 0.0f);
            mesh->v[idx].setZero();
            mesh->f[idx].setZero();
        }
    }

    // 2. 生成三角形索引
    for (int j = 0; j < resY; ++j) {
        for (int i = 0; i < resX; ++i) {
            uint32_t v0 = j * (resX + 1) + i;
            uint32_t v1 = v0 + 1;
            uint32_t v2 = (j + 1) * (resX + 1) + i;
            uint32_t v3 = v2 + 1;

            mesh->m_tris.emplace_back(v0, v1, v2);
            mesh->m_tris.emplace_back(v1, v3, v2);
        }
    }

    // 生成渲染用的表面索引 (对于 2D 布料，表面就是它自己)
    mesh->m_surface_tris = BuildSurfaceTriangles(mesh->m_tris);
}

void MeshBuilder::BuildBox(mesh_on_cpu* mesh, const float w, const float h, const float d, const int resX, const int resY, const int resZ) {
    const size_t num_nodes = (resX + 1) * (resY + 1) * (resZ + 1);
    PrepareMesh(mesh, num_nodes);

    const float dx = w / resX;
    const float dy = h / resY;
    const float dz = d / resZ;
    const Vec3 start = Vec3(-w/2, -h/2, -d/2);

    auto get_idx = [&](const int i, int j, int k) {
        return k * (resX + 1) * (resY + 1) + j * (resX + 1) + i;
    };

    // 1. 顶点
    for (int k = 0; k <= resZ; ++k) {
        for (int j = 0; j <= resY; ++j) {
            for (int i = 0; i <= resX; ++i) {
                mesh->p[get_idx(i,j,k)] = start + Vec3(i*dx, j*dy, k*dz);
            }
        }
    }

    // 2. 将每个小方块切分为 5 个四面体 (Tetrahedralization)
    for (int k = 0; k < resZ; ++k) {
        for (int j = 0; j < resY; ++j) {
            for (int i = 0; i < resX; ++i) {
                uint32_t v[8] = {
                    static_cast<uint32_t>(get_idx(i, j, k)),   static_cast<uint32_t>(get_idx(i + 1, j, k)),
                    static_cast<uint32_t>(get_idx(i + 1, j + 1, k)),   static_cast<uint32_t>(get_idx(i, j + 1, k)),
                    static_cast<uint32_t>(get_idx(i, j, k + 1)), static_cast<uint32_t>(get_idx(i + 1, j, k + 1)),
                    static_cast<uint32_t>(get_idx(i + 1, j + 1, k + 1)), static_cast<uint32_t>(get_idx(i, j + 1, k + 1))
                };
                // 经典 5-tet 分解法
                if ((i + j + k) % 2 == 0) {
                    mesh->m_tets.emplace_back(v[0], v[1], v[2], v[5]);
                    mesh->m_tets.emplace_back(v[0], v[2], v[3], v[7]);
                    mesh->m_tets.emplace_back(v[0], v[5], v[7], v[4]);
                    mesh->m_tets.emplace_back(v[2], v[5], v[7], v[6]);
                    mesh->m_tets.emplace_back(v[0], v[2], v[5], v[7]);
                } else {
                    mesh->m_tets.emplace_back(v[1], v[0], v[3], v[4]);
                    mesh->m_tets.emplace_back(v[1], v[3], v[2], v[6]);
                    mesh->m_tets.emplace_back(v[1], v[4], v[6], v[5]);
                    mesh->m_tets.emplace_back(v[3], v[4], v[6], v[7]);
                    mesh->m_tets.emplace_back(v[1], v[3], v[4], v[6]);
                }
            }
        }
    }
    // 自动提取外表面三角形供渲染
    mesh->m_surface_tris = BuildSurfaceTriangles(mesh->m_tets);
}

void MeshBuilder::BuildSphere(mesh_on_cpu* mesh, const float radius, const int sectors, const int stacks) {
    PrepareMesh(mesh, (sectors + 1) * (stacks + 1));

    // 1. 顶点生成 (UV 布局)
    for (int i = 0; i <= stacks; ++i) {
        const float phi = M_PI * static_cast<float>(i) / stacks;
        for (int j = 0; j <= sectors; ++j) {
            const float theta = 2.0f * M_PI * static_cast<float>(j) / sectors;
            float x = radius * sin(phi) * cos(theta);
            float y = radius * cos(phi);
            float z = radius * sin(phi) * sin(theta);

            mesh->p[i * (sectors + 1) + j] = Vec3(x, y, z);
        }
    }

    // 2. 索引生成
    for (int i = 0; i < stacks; ++i) {
        for (int j = 0; j < sectors; ++j) {
            uint32_t first = i * (sectors + 1) + j;
            uint32_t second = first + sectors + 1;

            mesh->m_tris.emplace_back(first, first + 1, second);
            mesh->m_tris.emplace_back(first + 1, second + 1, first);
        }
    }
    mesh->m_surface_tris = BuildSurfaceTriangles(mesh->m_tris);
}
