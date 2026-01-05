//
// Created by tianyan on 1/5/26.
//

#ifndef TAIYI_SHAPE_H
#define TAIYI_SHAPE_H

enum ShapeFlags {
    CollideShapes = 1u << 0,  // join shape-shape collide
    CollideParticles = 1u << 1,  // join particle-shape collide
    VisualOnly = 1u << 2,  // join no collide
};

struct Shape {

};

struct RigidBodyState {


};



#endif //TAIYI_SHAPE_H