//
// Created by tianyan on 12/23/25.
//

#include "MeshBuilder.h"

#ifdef _WIN64
float M_PI = 3.14159265358979323846;
#endif


void MeshBuilder::PrepareMesh(mesh_on_cpu* mesh, const size_t num_nodes) {
    if (!mesh) return;
    mesh->clear_topology(); // 清空旧的索引数据
    mesh->resize(num_nodes); // 调整 SoA 数组大小
    mesh->m_tris.reserve(num_nodes * 2);
    mesh->m_edges.reserve(num_nodes* 2);
    mesh->m_tets.reserve(num_nodes * 2);
}

/*
 * TODO: Add assert for number of node: should less than INVALID_VERTEX_ID
 */
void MeshBuilder::BuildCloth(mesh_on_cpu* mesh,
                             const float width, const float height,
                             const int resX, const int resY,
                             const Vec3& center,
                             const ClothOrientation orientation) {
    const size_t num_nodes = static_cast<size_t>(resX + 1) * static_cast<size_t>(resY + 1);
    PrepareMesh(mesh, num_nodes);

    const float dx = width  / static_cast<float>(resX);
    const float dy = height / static_cast<float>(resY);

    // 定义局部坐标系的基向量
    Vec3 u_dir, v_dir;
    if (orientation == ClothOrientation::Vertical) {
        u_dir = Vec3(1.0f, 0.0f, 0.0f);
        v_dir = Vec3(0.0f, 1.0f, 0.0f);
    } else {
        u_dir = Vec3(1.0f, 0.0f, 0.0f);
        v_dir = Vec3(0.0f, 0.0f, 1.0f);
    }

    const Vec3 start_pos = center - u_dir * (width * 0.5f) - v_dir * (height * 0.5f);

    // 1) 生成顶点
    for (int j = 0; j <= resY; ++j) {
        for (int i = 0; i <= resX; ++i) {
            const int idx = j * (resX + 1) + i;

            mesh->pos[idx]   = start_pos
                             + u_dir * (static_cast<float>(i) * dx)
                             + v_dir * (static_cast<float>(j) * dy);

            mesh->vel[idx].setZero();
            mesh->accel[idx].setZero();
        }
    }

    // 2) 生成三角形（索引 + 坐标）
    mesh->m_tris.clear();
    mesh->m_tris.reserve(static_cast<size_t>(resX) * static_cast<size_t>(resY) * 2);

    for (int j = 0; j < resY; ++j) {
        for (int i = 0; i < resX; ++i) {
            const auto v0 = static_cast<VertexID>( j      * (resX + 1) + i );
            const auto v1 = static_cast<VertexID>( v0 + 1 );
            const auto v2 = static_cast<VertexID>( (j + 1) * (resX + 1) + i );
            const auto v3 = static_cast<VertexID>( v2 + 1 );

            // 第一片： (v0, v1, v2)
            mesh->m_tris.emplace_back(
                v0, v1, v2,
                mesh->pos[v0], mesh->pos[v1], mesh->pos[v2]
            );

            // 第二片： (v1, v3, v2)
            mesh->m_tris.emplace_back(
                v1, v3, v2,
                mesh->pos[v1], mesh->pos[v3], mesh->pos[v2]
            );
        }
    }

    mesh->m_tris.shrink_to_fit();
    mesh->m_edges.shrink_to_fit();
    mesh->m_surface_tris = BuildSurfaceTriangles(mesh->m_tris);
}

void MeshBuilder::BuildBox(mesh_on_cpu* mesh, const float w, const float h, const float d) {
    ;
}

void MeshBuilder::BuildSphere(mesh_on_cpu* mesh, const float radius, const int sectors, const int stacks) {
    ;
}
