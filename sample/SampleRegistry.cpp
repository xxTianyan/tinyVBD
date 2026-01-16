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
    reg.Register(SampleId::EMPTY_SCENE, "Dummy Sample",
        []() -> SamplePtr { return std::make_unique<Sample>(); });

    reg.Register(SampleId::BASIC_CLOTH_EXAMPLE, "Basic Cloth Example",
        []()->SamplePtr { return std::make_unique<BasicCloth>();});

    reg.Register(SampleId::BUNNY, "Falling Bunny",
    []()->SamplePtr { return std::make_unique<FallingBunny>();});
}
