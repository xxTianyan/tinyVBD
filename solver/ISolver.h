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

    virtual void step(const float dt) = 0;
    virtual void init(const Scene& scene);
};


#endif //TAIYI_ISOLVER_H