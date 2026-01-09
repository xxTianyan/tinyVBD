//
// Created by tianyan on 1/9/26.
//

#ifndef TAIYI_BUNNY_H
#define TAIYI_BUNNY_H

#include "Sample.h"

class FallingBunny final : public Sample {
    void CreateWorld(AppContext &ctx) override;
    void Step(float dt) override;
    void BindShaders(AppContext &ctx) override;

private:
    size_t m_bunny_id_{};
};



#endif //TAIYI_BUNNY_H