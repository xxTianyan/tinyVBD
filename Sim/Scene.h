//
// Created by 徐天焱 on 2025/11/11.
//

#ifndef TAIYI_SCENE_H
#define TAIYI_SCENE_H

#include <span>

#include "Contacts.h"
#include "Builder.h"
#include "Model.h"
#include "Types.h"
struct MMaterial;

// all static values are stored in model struct, all changing values are stored in state struct

class Scene {
public:
    explicit Scene(MModel&& model): model_(std::move(model)) {
        state_in_ = model_.MakeState();
        state_out_ = model_.MakeState();
    };

    void InitStep();

    MModel model_;

    State state_in_;

    State state_out_;

    static bool RayNormal;

    void SetGravity(const Vec3& gravity){ gravity_ = gravity; };

private:
    Vec3 gravity_;
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
