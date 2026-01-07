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
size_t Builder::add_cloth(const float width, const float height,
                             const int resX, const int resY,
                             const Vec3& center,
                             const ClothOrientation orientation,
                             const char* name) const {

    if (resX <= 0 || resY <= 0)
        throw std::runtime_error("Builder::add_cloth: resX <= 0 || resY <= 0");

    const size_t particle_count = static_cast<size_t>(resX + 1) * static_cast<size_t>(resY + 1);

    // raylib Mesh indices is u16：the number of a single mesh must be less than 0xFFFF
    if (particle_count > static_cast<size_t>(std::numeric_limits<unsigned short>::max()) + 1ull) {
        throw std::runtime_error("Builder::add_cloth: particle_count > 65536, raylib u16 indices not supported. Split mesh.");
    }
    PrepareCapacity(particle_count);

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

    // add mesh info
    AddMeshInfo(name, particle_count, model_.edges.size(), model_.tris.size(), model_.tets.size());

    // assign mass and fix
    std::fill(model_.particle_inv_mass.begin(), model_.particle_inv_mass.end(), 1.0f);
    auto fix_left_z = [](const Vec3& pos) {
        if (std::abs(pos.z() - 1.0f) < 1e-5 ) return true;
        return false;
    };
    for (size_t i = model_.mesh_infos.back().particle.begin ; i < model_.mesh_infos.back().particle.end(); i++) {
        if (fix_left_z(model_.particle_pos0[i]))
            model_.particle_inv_mass[i] = 0.0f;
    }


    model_.topology_version++;

    return static_cast<int>(model_.mesh_infos.size()) - 1;
}


void Builder::PrepareCapacity(const size_t num) const {
    /*
     * TODO: Need to set different size of topology vec size in the future.
     */
    ensure_capacity(model_.particle_pos0, num);
    ensure_capacity(model_.particle_vel0, num);
    ensure_capacity(model_.particle_inv_mass, num);
    ensure_capacity(model_.edges, num*2);
    ensure_capacity(model_.tris, num*2);
    ensure_capacity(model_.tets, num*2);
    model_.num_particles += num;

    // resize std::vector<Vec3>, std::vector<float>
    model_.particle_pos0.resize(model_.num_particles);
    model_.particle_vel0.resize(model_.num_particles);
    model_.particle_inv_mass.resize(model_.num_particles);
}

void Builder::AddMeshInfo(const char* name, const size_t n_particle, const size_t n_edge,
            const size_t n_tri, const size_t n_tet) const {

    MeshInfo info{};
    info.name = name;
    info.particle = range{model_.particle_pos0.size() - n_particle, n_particle};
    info.edge = range{model_.edges.size() - n_edge, n_edge};
    info.tri = range{model_.tris.size() - n_tri, n_tri};
    info.tet = range{model_.tets.size() - n_tet, n_tet};

    model_.mesh_infos.emplace_back(info);
}


