//
// Created by xumiz on 2026/1/1.
//

#ifndef TINYVBD_BASIC_CLOTH_EXAMPLE_H
#define TINYVBD_BASIC_CLOTH_EXAMPLE_H
#include "Sample.h"

class HangingCloth final : public Sample {
public:
    HangingCloth() = default;
    void CreateWorld() override;
    void BindShaders() const;
};



#endif //TINYVBD_BASIC_CLOTH_EXAMPLE_H