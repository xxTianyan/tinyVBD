//
// Created by xumiz on 2026/1/2.
//

#include "Sample.h"
#include "SampleRegistry.h"


void RegisterAllSamples(SampleRegistry& reg) {

    // Empty Scene
    reg.Register(SampleId::DUMMY_SAMPLE, "Dummy Sample",
        []() -> SamplePtr { return std::make_unique<Sample>(); });

}
