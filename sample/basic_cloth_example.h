//
// Created by xumiz on 2026/1/1.
//

#ifndef TAIYI_BASIC_CLOTH_EXAMPLE_H
#define TAIYI_BASIC_CLOTH_EXAMPLE_H

#include "Sample.h"

class BasicCloth final : public Sample {
    void CreateWorld(AppContext &ctx) override;
    void Step(float dt) override;
    void BindShaders(AppContext &ctx) override;
    void Render(AppContext &ctx) override;

private:
    size_t m_cloth_id_{};
};




#endif //TAIYI_BASIC_CLOTH_EXAMPLE_H