//
// Created by xumiz on 2026/1/1.
//

#ifndef TINYVBD_SAMPLEREGISTRY_HPP
#define TINYVBD_SAMPLEREGISTRY_HPP

#include <functional>
#include <memory>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include "ISample.h"
#include "Types.h"

enum SampleId : int {
    empty_scene = 0,
    basic_cloth_sample = 1,
};

struct SampleInfo {
    SampleId id;
    const char* display_name;
};

class SampleRegistry {
public:
    using Factory = std::function<SamplePtr()>;

    void Register(const SampleId id, const char* display_name, Factory f) {
        infos_.push_back({id, display_name});
        factories_[id] = std::move(f);
    }

    SamplePtr Create(const SampleId id) const {
        const auto it = factories_.find(id);
        if (it == factories_.end()) throw std::runtime_error("Sample factory not found");
        return (it->second)();
    }

    const std::vector<SampleInfo>& Infos() const { return infos_; }

private:
    std::vector<SampleInfo> infos_;
    std::unordered_map<SampleId, Factory> factories_;
};


#endif //TINYVBD_SAMPLEREGISTRY_HPP