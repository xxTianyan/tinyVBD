//
// Created by xumiz on 2026/1/2.
//

#include "SampleRegistry.h"
#include "VBDDynamics.h"
#include "Sample.h"
#include "basic_cloth_example.h"

void RegisterAllSamples(SampleRegistry& reg) {

    // Empty Scene
    reg.Register(SampleId::DUMMY_SAMPLE, "Dummy Sample",
        []() -> SamplePtr { return std::make_unique<Sample>(); });

    reg.Register(SampleId::BASIC_CLOTH_EXAMPLE, "Basic Cloth Example",
        []()->SamplePtr { return std::make_unique<BasicCloth>();});

}
