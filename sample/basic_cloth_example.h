//
// Created by xumiz on 2026/1/1.
//

#ifndef BASIC_CLOTH_EXAMPLE_H
#define BASIC_CLOTH_EXAMPLE_H

#include "Sample.h"

class BasicCloth final : public Sample {
    void CreateWorld(AppContext &ctx) override;
    void Step(const float dt) override;
    void BindShaders(AppContext &ctx) override;
    void BuildRenderResources() override;
    void Render(AppContext &ctx) override;
};




#endif //BASIC_CLOTH_EXAMPLE_H