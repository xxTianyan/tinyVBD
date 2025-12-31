//
// Created by tianyan on 12/30/25.
//

#ifndef TINYVBD_MATERIALPARAMS_HPP
#define TINYVBD_MATERIALPARAMS_HPP

#include "Types.h"

enum MaterialType{
    StVK,
    NeoHookean,
};

struct MaterialID {
    MaterialType type;
    uint16_t id;
};

struct StVkMaterial {
    float E{};  // young's module
    float nu{};  // poisson's ratio
    float thickness = 1.0f;
};

inline StVkMaterial default_cloth() {
    return {1e6f, 0.3f};
};


#endif //TINYVBD_MATERIALPARAMS_HPP