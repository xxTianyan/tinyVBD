//
// Created by tianyan on 12/23/25.
//

#include "Builder.h"

#include "AdjacencyCSR.hpp"
#include "Model.h"

#ifdef _WIN64
float M_PI = 3.14159265358979323846;
#endif


size_t Builder::add_cloth(const float width, const float height,
                          const int resX, const int resY,
                          const Vec3& center,
                          const float mass,
                          const ClothOrientation orientation,
                          const FixSide fix_mask,
                          const char* name) const {
    // 1. Validation and basic setup
    if (resX <= 0 || resY <= 0)
        throw std::runtime_error("Builder::add_cloth: resX <= 0 || resY <= 0");
    if (width <= 0.0f || height <= 0.0f)
        throw std::runtime_error("Builder::add_cloth: width/height <= 0");

    // Total particles: (resX+1)*(resY+1) grid intersections + resX*resY quad centers
    const size_t grid_node_count = static_cast<size_t>(resX + 1) * static_cast<size_t>(resY + 1);
    const size_t center_node_count = static_cast<size_t>(resX) * static_cast<size_t>(resY);
    const size_t local_particle_count = grid_node_count + center_node_count;

    Builder::CheckVertexLimit(local_particle_count);

    const size_t base_particle = model_.mesh_infos.empty() ? 0ull : static_cast<size_t>(model_.mesh_infos.back().particle.end());
    PrepareCapacity(local_particle_count);

    const float dx = width  / static_cast<float>(resX);
    const float dy = height / static_cast<float>(resY);

    // 2. Setup local coordinate basis based on orientation
    Vec3 u_dir, v_dir;
    if (orientation == ClothOrientation::Vertical) {
        u_dir = Vec3(1.0f, 0.0f, 0.0f); v_dir = Vec3(0.0f, 1.0f, 0.0f);
    } else {
        u_dir = Vec3(1.0f, 0.0f, 0.0f); v_dir = Vec3(0.0f, 0.0f, 1.0f);
    }
    const Vec3 start_pos = center - u_dir * (width * 0.5f) - v_dir * (height * 0.5f);

    std::vector<float> mass_local(local_particle_count, 0.0f);

    // 3. Generate Grid Particles
    for (int j = 0; j <= resY; ++j) {
        for (int i = 0; i <= resX; ++i) {
            const size_t gid = base_particle + (j * (resX + 1) + i);
            model_.particle_pos0[gid] = start_pos + u_dir * (i * dx) + v_dir * (j * dy);
            model_.particle_vel0[gid].setZero();
        }
    }
    // 4. Generate Quad Center Particles
    for (int j = 0; j < resY; ++j) {
        for (int i = 0; i < resX; ++i) {
            const size_t gid = base_particle + grid_node_count + (j * resX + i);
            model_.particle_pos0[gid] = start_pos + u_dir * ((i + 0.5f) * dx) + v_dir * ((j + 0.5f) * dy);
            model_.particle_vel0[gid].setZero();
        }
    }

    const float total_mass = mass * static_cast<float>(local_particle_count);
    const float areal_density = total_mass / (width * height);

    auto tri_area = [](const Vec3& a, const Vec3& b, const Vec3& c) -> float {
        return 0.5f * (b - a).cross(c - a).norm();
    };

    // 5. Build Topology: Triangles and Bending Edges
    for (int j = 0; j < resY; ++j) {
        for (int i = 0; i < resX; ++i) {
            // Corner vertex indices
            const size_t i0 = j * (resX + 1) + i;
            const size_t i1 = i0 + 1;
            const size_t i2 = (j + 1) * (resX + 1) + i;
            const size_t i3 = i2 + 1;
            // Center vertex index
            const size_t iC = grid_node_count + (j * resX + i);

            const auto v0 = base_particle + i0;
            const auto v1 = base_particle + i1;
            const auto v2 = base_particle + i2;
            const auto v3 = base_particle + i3;
            const auto vC = base_particle + iC;

            const Vec3 &p0 = model_.particle_pos0[v0], &p1 = model_.particle_pos0[v1],
                       &p2 = model_.particle_pos0[v2], &p3 = model_.particle_pos0[v3],
                       &pC = model_.particle_pos0[vC];

            // Define 4 sub-triangles (counter-clockwise)
            struct TriDef { VertexID a, b, c; Vec3 pa, pb, pc; };
            TriDef subs[4] = {
                {v0, v1, vC, p0, p1, pC}, // Bottom Tri
                {v1, v3, vC, p1, p3, pC}, // Right Tri
                {v3, v2, vC, p3, p2, pC}, // Top Tri
                {v2, v0, vC, p2, p0, pC}  // Left Tri
            };

            for (const auto& t : subs) {
                model_.tris.emplace_back(t.a, t.b, t.c, t.pa, t.pb, t.pc);
                const float m_sub = areal_density * tri_area(t.pa, t.pb, t.pc);
                const float lumped_m = m_sub / 3.0f;
                mass_local[t.a - base_particle] += lumped_m;
                mass_local[t.b - base_particle] += lumped_m;
                mass_local[t.c - base_particle] += lumped_m;
            }

            // --- A) Internal Quad Bending (4 edges per quad) ---
            // These ensure the 4 internal triangles stay coplanar.
            // Edge v1-vC: Between Bottom and Right Tris
            model_.edges.emplace_back(v0, v3, v1, vC, p0, p3, p1, pC);
            // Edge v3-vC: Between Right and Top Tris
            model_.edges.emplace_back(v1, v2, v3, vC, p1, p2, p3, pC);
            // Edge v2-vC: Between Top and Left Tris
            model_.edges.emplace_back(v3, v0, v2, vC, p3, p0, p2, pC);
            // Edge v0-vC: Between Left and Bottom Tris
            model_.edges.emplace_back(v2, v1, v0, vC, p2, p1, p0, pC);

            // --- B) Inter-Quad Bending (Spanning across quads) ---
            // Horizontal neighbor: Shared edge v1-v3
            if (i < resX - 1) {
                const size_t iC_right = grid_node_count + (j * resX + (i + 1));
                const auto vC_right = static_cast<VertexID>(base_particle + iC_right);
                const Vec3& pC_right = model_.particle_pos0[vC_right];
                model_.edges.emplace_back(vC, vC_right, v1, v3, pC, pC_right, p1, p3);
            }
            // Vertical neighbor: Shared edge v2-v3
            if (j < resY - 1) {
                const size_t iC_up = grid_node_count + ((j + 1) * resX + i);
                const auto vC_up = static_cast<VertexID>(base_particle + iC_up);
                const Vec3& pC_up = model_.particle_pos0[vC_up];
                model_.edges.emplace_back(vC, vC_up, v2, v3, pC, pC_up, p2, p3);
            }
        }
    }

    // 6. Calculate counts for MeshInfo
    const size_t local_tri_count = 4 * resX * resY;
    const size_t internal_bending = 4 * resX * resY;
    const size_t horizontal_between = (resX - 1) * resY;
    const size_t vertical_between = resX * (resY - 1);
    const size_t local_edge_count = internal_bending + horizontal_between + vertical_between;

    AddMeshInfo(name, local_particle_count, local_edge_count, local_tri_count, 0);

    // 7. Assign mass and handle fixed boundaries
    for (size_t local_i = 0; local_i < local_particle_count; ++local_i) {
        const size_t gid = base_particle + local_i;
        const float m = mass_local[local_i];
        model_.particle_inv_mass[gid] = (m > 1e-12f) ? (1.0f / m) : 0.0f;

        // Only grid nodes (intersections) can be on the physical boundary
        if (local_i < grid_node_count) {
            const int i = static_cast<int>(local_i % (resX + 1));
            const int j = static_cast<int>(local_i / (resX + 1));

            bool should_fix = false;
            if ((fix_mask & TOP)    && (j == resY)) should_fix = true;
            if ((fix_mask & BOTTOM) && (j == 0))    should_fix = true;
            if ((fix_mask & LEFT)   && (i == 0))    should_fix = true;
            if ((fix_mask & RIGHT)  && (i == resX)) should_fix = true;

            if (should_fix) model_.particle_inv_mass[gid] = 0.0f;
        }
    }

    model_.topology_version++;
    return model_.mesh_infos.size() - 1;
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




