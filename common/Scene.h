//
// Created by 徐天焱 on 2025/11/11.
//

#ifndef TAIYI_SCENE_H
#define TAIYI_SCENE_H

#include <span>

#include "Contacts.h"
#include "Model.h"
#include "Types.h"
struct MMaterial;

struct SimView {
    std::span<Vec3> pos;
    std::span<Vec3> prev_pos;
    std::span<Vec3> inertia_pos;
    std::span<Vec3> vel;
    std::span<Vec3> accel;
    std::span<float> inv_mass;
    const std::span<uint8_t> fixed;
    const std::span<edge> edges;
    const std::span<triangle> tris;
    const std::span<tetrahedron> tets;
    const Vec3& gravity;
    const MMaterial& material_params;
};


// all static values are stored in model struct, all changing values are stored in state struct

struct SceneModel {
    Vec3 gravity;
    std::vector<MMaterial> materials;
    std::vector<MaterialID> mesh_to_material;

    std::vector<MModel> meshes;

    // ShapeStore shapes;
    size_t total_vertices;
};

struct SceneState {
    std::vector<State> meshes;

    void BeginStep();
};

class Scene {
public:
    Scene() = default;

    void InitStep();

    MeshID Add(MModel&& model);

    SimView MakeSimView(MeshID id);

    //
    [[nodiscard]] const MeshModels& GetMeshModels() const { return models; }
    States states_in;

    // temporary
    static bool RayNormal;

private:
    Vec3 gravity;
    MeshModels models;
};

// helper functions for constructing scene
template <class Vec>
size_t checked_index(const int32_t id, const Vec& v, const char* what) {
    if (id < 0) {
        throw std::runtime_error(std::string(what) + " id < 0");
    }
    const auto idx = static_cast<size_t>(id);
    if (idx >= v.size()) {
        throw std::runtime_error(std::string(what) + " id out of range");
    }
    return idx;
}




























#endif //TAIYI_SCENE_H
