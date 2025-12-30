//
// Created by tianyan on 12/30/25.
//

#ifndef TINYVBD_MATERIALPARAMS_HPP
#define TINYVBD_MATERIALPARAMS_HPP

#include "Types.h"

struct StVkMaterial {
    float lambda; // first lame parameter
    float mu;  // second lame parameter
    float thickness;
};


#endif //TINYVBD_MATERIALPARAMS_HPP