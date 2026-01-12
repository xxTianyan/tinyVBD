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

enum FixSide { NONE = 0, TOP = 1, BOTTOM = 2, LEFT = 4, RIGHT = 8 };

class Builder {

public:
    explicit Builder(MModel& model) : model_(model) {};

    [[nodiscard]] size_t add_cloth(float width, float height, int resX, int resY, const Vec3& center = Vec3{0.0f,0.0f,0.0f},
                        float mass = .1f, ClothOrientation orientation = ClothOrientation::Horizontal, FixSide fix_mask = FixSide::TOP, const char* = "cloth") const;
    [[nodiscard]] size_t add_bunny(float mass = 1.f) const;

    [[nodiscard]] size_t add_single_tet() const;

    // void add_rigidbody();

private:

    MModel& model_;

    void PrepareCapacity(size_t num) const;
    // void PrepareTopologyCapacity(size_t num) const;

    void AddMeshInfo(const char* name, size_t n_particle, size_t n_edge,
                size_t n_tri, size_t n_tet) const;

    static void CheckVertexLimit(const uint32_t local_particle_count) {
        if (local_particle_count > static_cast<size_t>(std::numeric_limits<unsigned short>::max()) + 1ull) {
            throw std::runtime_error("Builder: particle_count > 65536, raylib u16 indices not supported.");
        }
    }
};



#endif //TAIYI_BUILDER_H