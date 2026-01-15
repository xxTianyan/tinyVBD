//
// Created by tianyan on 12/29/25.
//

#ifndef TAIYI_ADJACENCYCSR_HPP
#define TAIYI_ADJACENCYCSR_HPP

#include <vector>
#include "Types.h"

struct AdjacencyCSR {
    // offsets.size() = num_nodes + 1
    std::vector<uint32_t> offsets;

    // packed incident entries, length == offsets.back()
    // packed = (element_id << 2) | local_order
    std::vector<uint32_t> incidents;

    static constexpr uint32_t ORDER_MASK = 0x3u;  // low 2 bits
    static inline uint32_t pack(const uint32_t element_id, const uint32_t local_order) {
        return (element_id << 2) | (local_order & ORDER_MASK);
    }
    static inline uint32_t unpack_id(const uint32_t packed) {
        return packed >> 2;
    }
    static inline uint32_t unpack_order(const uint32_t packed) {
        return packed & ORDER_MASK;
    }

    [[nodiscard]] inline uint32_t begin(const VertexID v) const { return offsets[v]; }
    [[nodiscard]] inline uint32_t end(const VertexID v) const { return offsets[v + 1]; }
    [[nodiscard]] inline uint32_t degree(const VertexID v) const { return end(v) - begin(v); }
    [[nodiscard]] inline uint32_t incident_packed(const VertexID v, const uint32_t k) const { return incidents[begin(v) + k]; }

};

struct ForceElementAdjacencyInfo {
    AdjacencyCSR vertex_edges;
    AdjacencyCSR vertex_faces;
    AdjacencyCSR vertex_tets;
};

template <class Elem, class GetVertex>
void BuildVertexIncidentCSR(const size_t num_nodes,
    const std::vector<Elem>& elems,
    const uint32_t verts_per_elem,
    GetVertex get_vertex,
    AdjacencyCSR& adj) {
    auto& offsets = adj.offsets;
    auto& incidents = adj.incidents;

    offsets.assign(num_nodes + 1, 0u);
    incidents.clear();

    if (elems.empty()) {
        return;
    }

    for (uint32_t elem_id = 0; elem_id < static_cast<uint32_t>(elems.size()); ++elem_id) {
        const auto& elem = elems[elem_id];
        for (uint32_t k = 0; k < verts_per_elem; ++k) {
            const auto v = static_cast<uint32_t>(get_vertex(elem, k));
            offsets[v + 1] += 1u;
        }
    }

    for (size_t i = 1; i < offsets.size(); ++i) {
        offsets[i] += offsets[i - 1];
    }

    incidents.resize(offsets.back());
    std::vector<uint32_t> cursor = offsets;

    for (uint32_t elem_id = 0; elem_id < static_cast<uint32_t>(elems.size()); ++elem_id) {
        const auto& elem = elems[elem_id];
        for (uint32_t k = 0; k < verts_per_elem; ++k) {
            const auto v = static_cast<uint32_t>(get_vertex(elem, k));
            const uint32_t dst = cursor[v]++;
            incidents[dst] = AdjacencyCSR::pack(elem_id, k);
        }
    }
}


#endif //TAIYI_ADJACENCYCSR_HPP