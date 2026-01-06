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
                        ClothOrientation orientation = ClothOrientation::Vertical) const;

    // void add_rigidbody();

private:

    MModel& model_;

    void PrepareCapacity(const size_t num) const;


    void AddMeshInfo(const char* name, const size_t n_particle, const size_t n_edge,
                const size_t n_tri, const size_t n_tet) const;

};



#endif //TAIYI_BUILDER_H