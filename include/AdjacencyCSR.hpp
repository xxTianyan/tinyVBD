//
// Created by tianyan on 12/29/25.
//

#ifndef TINYVBD_ADJACENCYCSR_HPP
#define TINYVBD_ADJACENCYCSR_HPP

#include <vector>

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

    inline uint32_t begin(const VertexId v) const { return offsets[v]; }
    inline uint32_t end(const VertexId v) const { return offsets[v + 1]; }
    inline uint32_t degree(const VertexId v) const { return end(v) - begin(v); }
    inline uint32_t incident_packed(const VertexId v, const uint32_t k) const { return incidents[begin(v) + k]; }

};

struct ForceElementAdjacencyInfo {
    AdjacencyCSR vertex_edges;
    AdjacencyCSR vertex_faces;
    AdjacencyCSR vertex_tets;
};




#endif //TINYVBD_ADJACENCYCSR_HPP