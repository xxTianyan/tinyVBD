//
// Created by 徐天焱 on 2025/11/11.
//

#ifndef TAIYI_SCENE_H
#define TAIYI_SCENE_H

#include <span>

#include "Contacts.h"
#include "Mesh.h"
#include "Types.h"
struct MMaterial;

struct SimView {
    std::span<Vec3> pos;
    std::span<Vec3> prev_pos;
    std::span<Vec3> inertia_pos;
    std::span<Vec3> vel;
    std::span<Vec3> accel;
    // std::span<Vec3> normal;
    std::span<float> inv_mass;
    const std::span<uint8_t> fixed;
    const std::span<edge> edges;
    const std::span<triangle> tris;
    const std::span<tetrahedron> tets;
    const Vec3& gravity;
    const MMaterial& material_params;
    const ForceElementAdjacencyInfo& adj;
};


// all static values are stored in model struct, all changing values are stored in state struct

struct SceneModel {
    Vec3 gravity;
    std::vector<MMaterial> materials;
    std::vector<MaterialID> mesh_to_material;

    std::vector<MeshModel> meshes;

    // ShapeStore shapes;
    size_t total_vertices;
};

struct SceneState {
    std::vector<MeshState> meshes;

    void BeginStep();
};

class Scene {
public:
    explicit Scene(const Vec3 &gravity);

    MeshID Add(MeshModel&& model, MeshState&& state);
    MaterialID AddMaterial(const MMaterial& material);

    SimView MakeSimView(MeshID mesh_id);

    // maybe need dt here
    void InitStep();

    // give up material part for now, since it may be moved to shape in the future.
    void BindMeshMaterial(MeshID mesh, MaterialID mat);

    // if some vertex need tobe fixed
    void ApplyFixConsition(MeshID mesh, const std::function<bool(const Vec3&)>& pred);

    // for rendering
    [[nodiscard]] const std::vector<MeshModel>& MeshModels() const { return model_.meshes; }
    [[nodiscard]] const std::vector<MeshState>& MeshStates() const { return state_.meshes; }
    std::vector<MeshState>& MeshStates() { return state_.meshes; }

    Contacts& ContactsRef() { return contacts_; }
    [[nodiscard]] const SceneModel& Model() const { return model_; }
    [[nodiscard]] const SceneState& State() const { return state_; }

    // temporary
    static bool RayNormal;

private:
    SceneModel model_;
    SceneState state_;
    Contacts contacts_;

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

static void validate_mesh_sizes(const MeshModel& mm, const MeshState& ms) {
    const size_t n = mm.size();

    auto req_eq = [&](const size_t sz, const char* name) {
        if (sz != n) {
            throw std::runtime_error(std::string("MeshState size mismatch: ") + name);
        }
    };

    req_eq(ms.prev_pos.size(),    "prev_pos");
    req_eq(ms.inertia_pos.size(), "inertia_pos");
    req_eq(ms.vel.size(),         "vel");
    req_eq(ms.accel.size(),       "accel");
    // req_eq(ms.n.size(),           "n");

    // Model-side per-vertex arrays (at minimum inv_mass/fixed)
    if (mm.inv_mass.size() != n) {
        throw std::runtime_error("MeshModel size mismatch: inv_mass");
    }
    if (mm.fixed.size() != n) {
        throw std::runtime_error("MeshModel size mismatch: fixed");
    }
}


























#endif //TAIYI_SCENE_H
