//
// Created by tianyan on 12/29/25.
//

#ifndef TAIYI_ADJACENCYCSR_HPP
#define TAIYI_ADJACENCYCSR_HPP

#include <vector>
#include <cstdint>
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

// adj test
/*const auto& mesh = falling_cloth.m_world->meshes[0];
const auto& tris = mesh->m_tris;
const auto& adj_info = mesh->adjacencyInfo;
const VertexId v = 177;
for (const auto& t: tris) {
    auto v1 = ToRayVec(mesh->pos[t.vertices[0]]);
    auto v2 = ToRayVec(mesh->pos[t.vertices[1]]);
    auto v3 = ToRayVec(mesh->pos[t.vertices[2]]);
    DrawLine3D(v1, v2, MAROON);
    DrawLine3D(v2, v3, MAROON);
    DrawLine3D(v3, v1, MAROON);
}
DrawSphere(ToRayVec(mesh->pos[v]), 0.01, BLUE);
for (uint32_t f = adj_info.vertex_faces.begin(v); f < adj_info.vertex_faces.end(v); ++f) {
    auto pack_id = adj_info.vertex_faces.incidents[f];
    auto face_id = AdjacencyCSR::unpack_id(pack_id);
    auto the_tri = mesh->m_tris[face_id];
    auto v1 = ToRayVec(mesh->pos[the_tri.vertices[0]]);
    auto v2 = ToRayVec(mesh->pos[the_tri.vertices[1]]);
    auto v3 = ToRayVec(mesh->pos[the_tri.vertices[2]]);
    DrawTriangle3D(v3, v2, v1, RAYWHITE);
}*/


#endif //TAIYI_ADJACENCYCSR_HPP