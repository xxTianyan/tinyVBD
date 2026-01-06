//
// Created by tianyan on 12/23/25.
//

#include "Builder.h"
#include "Model.h"

#ifdef _WIN64
float M_PI = 3.14159265358979323846;
#endif


/*
 * TODO: Add assert for number of node: should less than INVALID_VERTEX_ID
 */
MeshID Builder::add_cloth(const float width, const float height,
                             const int resX, const int resY,
                             const Vec3& center,
                             const ClothOrientation orientation) {

    if (resX <= 0 || resY <= 0)
        throw std::runtime_error("Builder::add_cloth: resX <= 0 || resY <= 0");

    const size_t num_nodes = static_cast<size_t>(resX + 1) * static_cast<size_t>(resY + 1);
    ResizeDeformable(num_nodes);

    const float dx = width  / static_cast<float>(resX);
    const float dy = height / static_cast<float>(resY);

    // local coordinates and basis vector
    Vec3 u_dir, v_dir;
    if (orientation == ClothOrientation::Vertical) {
        u_dir = Vec3(1.0f, 0.0f, 0.0f);
        v_dir = Vec3(0.0f, 1.0f, 0.0f);
    } else {
        u_dir = Vec3(1.0f, 0.0f, 0.0f);
        v_dir = Vec3(0.0f, 0.0f, 1.0f);
    }

    const Vec3 start_pos = center - u_dir * (width * 0.5f) - v_dir * (height * 0.5f);

    // 1) vertex
    for (int j = 0; j <= resY; ++j) {
        for (int i = 0; i <= resX; ++i) {
            const int idx = j * (resX + 1) + i;

            model_.particle_pos0[idx]   = start_pos
                             + u_dir * (static_cast<float>(i) * dx)
                             + v_dir * (static_cast<float>(j) * dy);

            model_.particle_vel0[idx].setZero();
        }
    }

    // 2) triangles and edges
    ReserveTopology(num_nodes * 2);

    for (int j = 0; j < resY; ++j) {
        for (int i = 0; i < resX; ++i) {
            const auto v0 = static_cast<VertexID>( j * (resX + 1) + i );
            const auto v1 = static_cast<VertexID>( v0 + 1 );
            const auto v2 = static_cast<VertexID>( (j + 1) * (resX + 1) + i );
            const auto v3 = static_cast<VertexID>( v2 + 1 );

            // (v0, v1, v2)
            model_.tris.emplace_back(
                v0, v1, v2,
                model_.particle_pos0[v0], model_.particle_pos0[v1], model_.particle_pos0[v2]
            );

            // (v1, v3, v2)
            model_.tris.emplace_back(
                v1, v3, v2,
                model_.particle_pos0[v1], model_.particle_pos0[v3], model_.particle_pos0[v2]
            );

            // (v0, v3, v1, v2)
            model_.edges.emplace_back(
                v0, v3, v1, v2,
                model_.particle_pos0[v0], model_.particle_pos0[v3], model_.particle_pos0[v1], model_.particle_pos0[v2]);
        }
    }

    ShrinkTopology();

    // add mesh info


    return 1;
}




