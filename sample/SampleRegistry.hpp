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
#include "Sample.h"


enum SampleId : int {
    basic_cloth_sample = 0,
};

struct SampleInfo {
    SampleId id;
    const char* display_name;
};

class SampleRegistry {
public:
    using Factory = std::function<std::unique_ptr<ISample>()>;

    void Register(const SampleId id, const char* display_name, Factory f) {
        infos_.push_back({id, display_name});
        factories_[id] = std::move(f);
    }

    std::unique_ptr<ISample> Create(const SampleId id) const {
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