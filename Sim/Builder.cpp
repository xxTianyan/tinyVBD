//
// Created by tianyan on 12/23/25.
//

#include "Builder.h"
#include "Model.h"

#ifdef _WIN64
float M_PI = 3.14159265358979323846;
#endif


size_t Builder::add_cloth(const float width, const float height,
                             const int resX, const int resY,
                             const Vec3& center,
                             const float density,
                             const ClothOrientation orientation,
                             const char* name) const {
    {
    if (resX <= 0 || resY <= 0)
        throw std::runtime_error("Builder::add_cloth: resX <= 0 || resY <= 0");
    if (width <= 0.0f || height <= 0.0f)
        throw std::runtime_error("Builder::add_cloth: width/height <= 0");
    if (density <= 0.0f)
        throw std::runtime_error("Builder::add_cloth: density <= 0 (expect kg/m^2)");

    const size_t local_particle_count =
        static_cast<size_t>(resX + 1) * static_cast<size_t>(resY + 1);

    // raylib Mesh indices is u16: single mesh must be < 65536 vertices (indexable by u16)
    if (local_particle_count > static_cast<size_t>(std::numeric_limits<unsigned short>::max()) + 1ull) {
        throw std::runtime_error("Builder::add_cloth: particle_count > 65536, raylib u16 indices not supported. Split mesh.");
    }

    // ----- base offsets -----
    const size_t base_particle =
        model_.mesh_infos.empty() ? 0ull : static_cast<size_t>(model_.mesh_infos.back().particle.end());


    PrepareCapacity(local_particle_count);

    const float dx = width  / static_cast<float>(resX);
    const float dy = height / static_cast<float>(resY);

    // local basis
    Vec3 u_dir, v_dir;
    if (orientation == ClothOrientation::Vertical) {
        u_dir = Vec3(1.0f, 0.0f, 0.0f);
        v_dir = Vec3(0.0f, 1.0f, 0.0f);
    } else {
        u_dir = Vec3(1.0f, 0.0f, 0.0f);
        v_dir = Vec3(0.0f, 0.0f, 1.0f);
    }

    const Vec3 start_pos = center - u_dir * (width * 0.5f) - v_dir * (height * 0.5f);

    // per-vertex mass accumulator (local indexing)
    std::vector<float> mass_local(local_particle_count, 0.0f);

    // 1) vertices
    for (int j = 0; j <= resY; ++j) {
        for (int i = 0; i <= resX; ++i) {
            const size_t local_idx = static_cast<size_t>(j) * static_cast<size_t>(resX + 1) + static_cast<size_t>(i);
            const size_t gid = base_particle + local_idx;

            model_.particle_pos0[gid] =
                start_pos
                + u_dir * (static_cast<float>(i) * dx)
                + v_dir * (static_cast<float>(j) * dy);

            model_.particle_vel0[gid].setZero();
        }
    }

    // 2) triangles and edges (+ accumulate lumped mass from triangle areas)
    auto tri_area = [](const Vec3& a, const Vec3& b, const Vec3& c) -> float {
        return 0.5f * (b - a).cross(c - a).norm();
    };

    for (int j = 0; j < resY; ++j) {
        for (int i = 0; i < resX; ++i) {
            const auto i0 = static_cast<uint32_t>(j * (resX + 1) + i);
            const uint32_t i1 = i0 + 1;
            const auto i2 = static_cast<uint32_t>((j + 1) * (resX + 1) + i);
            const uint32_t i3 = i2 + 1;

            const auto v0 = static_cast<VertexID>(base_particle + i0);
            const auto v1 = static_cast<VertexID>(base_particle + i1);
            const auto v2 = static_cast<VertexID>(base_particle + i2);
            const auto v3 = static_cast<VertexID>(base_particle + i3);

            const Vec3& p0 = model_.particle_pos0[v0];
            const Vec3& p1 = model_.particle_pos0[v1];
            const Vec3& p2 = model_.particle_pos0[v2];
            const Vec3& p3 = model_.particle_pos0[v3];

            // (v0, v1, v2)
            model_.tris.emplace_back(v0, v1, v2, p0, p1, p2);
            {
                const float A = tri_area(p0, p1, p2);
                const float m = density * A;
                const float add = m / 3.0f;
                mass_local[i0] += add;
                mass_local[i1] += add;
                mass_local[i2] += add;
            }

            // (v1, v3, v2)
            model_.tris.emplace_back(v1, v3, v2, p1, p3, p2);
            {
                const float A = tri_area(p1, p3, p2);
                const float m = density * A;
                const float add = m / 3.0f;
                mass_local[i1] += add;
                mass_local[i3] += add;
                mass_local[i2] += add;
            }

            // bending edge (v0, v3, v1, v2)  -- keep your original convention
            model_.edges.emplace_back(v0, v3, v1, v2, p0, p3, p1, p2);
        }
    }

    // add mesh info (kept as your original call pattern)
    AddMeshInfo(name, local_particle_count, model_.edges.size(), model_.tris.size(), model_.tets.size());

    // 3) assign inverse mass from lumped masses (only for this cloth's vertex range)
    //    inv_mass = 1 / mass (kg^-1). fixed => inv_mass = 0.

    // Example: fix the "top edge" along +v_dir (invariant to center/orientation)
    auto is_fixed = [&](const Vec3& p) -> bool {
        // coordinate along v_dir relative to center should be +height/2
        const float s = (p - center).dot(v_dir); // v_dir is unit
        return std::abs(s - height * 0.5f) < 1e-5f;
    };

    for (size_t local_i = 0; local_i < local_particle_count; ++local_i) {
        constexpr float eps_mass = 1e-12f;
        const size_t gid = base_particle + local_i;

        const float m = mass_local[local_i];
        model_.particle_inv_mass[gid] = (m > eps_mass) ? (1.0f / m) : 0.0f;

        if (is_fixed(model_.particle_pos0[gid])) {
            model_.particle_inv_mass[gid] = 0.0f;
        }
    }

    model_.topology_version++;

    return model_.mesh_infos.size() - 1;
}
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


