//
// Created by xumiz on 2026/1/2.
//

#include "SampleRegistry.h"
#include "VBDSolver.h"
#include "Sample.h"
#include "hanging_cloth.hpp"
#include "falling_bunny.hpp"

void RegisterAllSamples(SampleRegistry& reg) {

    // Empty Scene
    reg.Register(SampleId::EMPTY_SCENE, "Empty Scene",
        []() -> SamplePtr { return std::make_unique<Sample>(); });

    reg.Register(SampleId::HANGING_CLOTH, "Hanging Cloth",
        []()->SamplePtr { return std::make_unique<HangingCloth>();});

    reg.Register(SampleId::FALLING_BUNNY, "Falling Bunny",
    []()->SamplePtr { return std::make_unique<FallingBunny>();});
}
