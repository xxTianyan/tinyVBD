//
// Created by tianyan on 1/6/26.
//

#ifndef TAIYI_ISOLVER_H
#define TAIYI_ISOLVER_H

class Scene;

class ISolver {
public:
    ISolver() = default;
    virtual ~ISolver() = default;

    virtual void Init(const Scene& scene) = 0;
    virtual void Step(Scene& scene, float dt) = 0;
};


#endif //TAIYI_ISOLVER_H
