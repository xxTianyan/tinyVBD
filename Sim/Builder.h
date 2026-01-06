//
// Created by tianyan on 12/23/25.
//

#ifndef TAIYI_BUILDER_H
#define TAIYI_BUILDER_H

#include "Model.h"
#include "raylib.h"
#include "Types.h"

struct MModel;
struct mesh_on_cpu;

template<class T>
static void ensure_capacity(std::vector<T>& v, size_t extra, const size_t min_grow = 1024) {
    const size_t need = v.size() + extra;
    if (need <= v.capacity()) return;

    size_t new_cap = v.capacity();
    if (new_cap == 0) new_cap = min_grow;

    while (new_cap < need) {
        const size_t prev = new_cap;
        new_cap = new_cap + (new_cap >> 1); // 1.5x
        if (new_cap <= prev) {
            new_cap = need;
            break;
        }
    }
    v.reserve(new_cap);
}



enum class ClothOrientation {
    Vertical,   //  (XY Plane)
    Horizontal  //  (XZ Plane)
};

class Builder {

public:
    explicit Builder(MModel& model) : model_(model) {};

    MeshID add_cloth(float width, float height, int resX, int resY, const Vec3& center = Vec3(0,0,0),
                        ClothOrientation orientation = ClothOrientation::Vertical);

    // void add_rigidbody();

private:

    MModel& model_;

    void PrepareCapcity(const size_t num) {
        ensure_capacity(model_.particle_pos0, num);
        ensure_capacity(model_.particle_vel0, num);
        ensure_capacity(model_.particle_inv_mass, num);
        ensure_capacity(model_.edges, num*2);
        ensure_capacity(model_.tris, num*2);
        ensure_capacity(model_.tets, num*2);
        model_.num_particles += num;
    }

    void ResizeDeformable(const size_t n_nodes) const {
        model_.num_particles += n_nodes;
        const auto num = model_.num_particles;
        model_.particle_inv_mass.resize(num);
        model_.if_particle_fixed.resize(num);
        model_.particle_pos0.resize(num);
        model_.particle_vel0.resize(num);
    }
    void ReserveTopology(const size_t num) const {
        model_.edges.reserve(num);
        model_.tris.reserve(num);
        model_.tets.reserve(num);
    }


    void AddMeshInfo(const char* name, const size_t n_particle, const size_t n_edge,
                const size_t n_tri, const size_t n_tet) const {

        MeshInfo info{};
        model_.name_pool_.emplace_back(name);
        info.name = model_.name_pool_.back().c_str();

        info.particle = range{model_.particle_pos0.size(), n_particle};
        info.edge = range{model_.edges.size(), n_edge};
        info.tri = range{model_.tris.size(), n_tri};
        info.tet = range{model_.tets.size(), n_tet};

    }

};



#endif //TAIYI_BUILDER_H