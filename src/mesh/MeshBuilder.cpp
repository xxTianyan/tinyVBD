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
}

void MeshBuilder::BuildCloth(mesh_on_cpu* mesh,
                             const float width, const float height,
                             const int resX, const int resY,
                             const Vec3& center,
                             const ClothOrientation orientation)
{
    const size_t num_nodes = (resX + 1) * (resY + 1);
    PrepareMesh(mesh, num_nodes);

    const float dx = width / resX;
    const float dy = height / resY;


    // 定义局部坐标系的基向量
    Vec3 u_dir, v_dir;

    if (orientation == ClothOrientation::Vertical) {
        // 竖直模式 (XY平面): 宽沿X轴，高沿Y轴
        u_dir = Vec3(1.0f, 0.0f, 0.0f);
        v_dir = Vec3(0.0f, 1.0f, 0.0f);
    } else {
        // 水平模式 (XZ平面): 宽沿X轴，高沿Z轴
        u_dir = Vec3(1.0f, 0.0f, 0.0f);
        v_dir = Vec3(0.0f, 0.0f, 1.0f);
    }

    // 计算起始点 (左下角/左上角，取决于你的纹理坐标习惯，这里保持几何中心对齐)
    // start_pos = center - (宽的一半 * u方向) - (高的一半 * v方向)
    const Vec3 start_pos = center - u_dir * (width * 0.5f) - v_dir * (height * 0.5f);

    // 1. 生成顶点
    for (int j = 0; j <= resY; ++j) {
        for (int i = 0; i <= resX; ++i) {
            const int idx = j * (resX + 1) + i;

            // 现在的顶点位置 = 起点 + (i * dx * u方向) + (j * dy * v方向)
            mesh->pos[idx] = start_pos + u_dir * (float(i) * dx) + v_dir * (float(j) * dy);

            mesh->vel[idx].setZero();
            mesh->accel[idx].setZero();

            // 如果你的 mesh 结构体里有法线(normal)，别忘了在这里初始化
            // Vertical 模式法线通常是 (0, 0, 1) 或 (0, 0, -1)
            // Horizontal 模式法线通常是 (0, 1, 0)
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

    mesh->m_surface_tris = BuildSurfaceTriangles(mesh->m_tris);
}

// MeshBuilder.h 里建议新增这个重载
// static void BuildBox(mesh_on_cpu* mesh, float w, float h, float d);

inline float SignedTetVolume6(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    // 6 * volume (signed)
    return (b - a).dot((c - a).cross(d - a));
}

void MeshBuilder::BuildBox(mesh_on_cpu* mesh, const float w, const float h, const float d) {
    // 0) 基本检查
    if (!mesh) return;
    if (w <= 0.0f || h <= 0.0f || d <= 0.0f) return;

    // 1) 顶点：8 个角点
    PrepareMesh(mesh, 8);

    // 如果 PrepareMesh 不会清拓扑，这里显式清一下更安全
    mesh->m_tets.clear();
    mesh->m_surface_tris.clear();

    const Vec3 p0(-w * 0.5f, -h * 0.5f, -d * 0.5f);
    const Vec3 p1(+w * 0.5f, -h * 0.5f, -d * 0.5f);
    const Vec3 p2(+w * 0.5f, +h * 0.5f, -d * 0.5f);
    const Vec3 p3(-w * 0.5f, +h * 0.5f, -d * 0.5f);
    const Vec3 p4(-w * 0.5f, -h * 0.5f, +d * 0.5f);
    const Vec3 p5(+w * 0.5f, -h * 0.5f, +d * 0.5f);
    const Vec3 p6(+w * 0.5f, +h * 0.5f, +d * 0.5f);
    const Vec3 p7(-w * 0.5f, +h * 0.5f, +d * 0.5f);

    mesh->pos[0] = p0; mesh->pos[1] = p1; mesh->pos[2] = p2; mesh->pos[3] = p3;
    mesh->pos[4] = p4; mesh->pos[5] = p5; mesh->pos[6] = p6; mesh->pos[7] = p7;

    // 2) 体积剖分：标准 5-tet（沿体对角线 0-6）
    auto add_tet_ccw = [&](uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
        // 保证四面体有正体积（一致取向），避免后续算法/抽面/法线混乱
        const Vec3& A = mesh->pos[a];
        const Vec3& B = mesh->pos[b];
        const Vec3& C = mesh->pos[c];
        const Vec3& D = mesh->pos[d];
        float v6 = SignedTetVolume6(A, B, C, D);
        if (v6 < 0.0f) std::swap(c, d);
        mesh->m_tets.emplace_back(a, b, c, d);
    };

    add_tet_ccw(0, 1, 2, 6);
    add_tet_ccw(0, 2, 3, 6);
    add_tet_ccw(0, 3, 7, 6);
    add_tet_ccw(0, 7, 4, 6);
    add_tet_ccw(0, 4, 5, 6);

    // 3) 渲染外表面：直接生成 6 个面 * 2 三角形，保证 outward winding
    auto add_tri = [&](uint32_t a, uint32_t b, uint32_t c) {
        mesh->m_surface_tris.push_back(a);
        mesh->m_surface_tris.push_back(b);
        mesh->m_surface_tris.push_back(c);
    };

    // -Z 面（底面），外法线 -Z：四边形 [0,1,2,3] -> 用 (0,3,2) (0,2,1)
    add_tri(0, 3, 2);
    add_tri(0, 2, 1);

    // +Z 面（顶面），外法线 +Z： [4,5,6,7] -> (4,5,6) (4,6,7)
    add_tri(4, 5, 6);
    add_tri(4, 6, 7);

    // -Y 面，外法线 -Y： [0,1,5,4] -> (0,1,5) (0,5,4)
    add_tri(0, 1, 5);
    add_tri(0, 5, 4);

    // +Y 面，外法线 +Y： [3,2,6,7] -> (3,7,6) (3,6,2)
    add_tri(3, 7, 6);
    add_tri(3, 6, 2);

    // -X 面，外法线 -X： [0,3,7,4] -> (0,4,7) (0,7,3)
    add_tri(0, 4, 7);
    add_tri(0, 7, 3);

    // +X 面，外法线 +X： [1,2,6,5] -> (1,2,6) (1,6,5)
    add_tri(1, 2, 6);
    add_tri(1, 6, 5);
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

            mesh->pos[i * (sectors + 1) + j] = Vec3(x, y, z);
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
