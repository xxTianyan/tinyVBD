//
// Created by tianyan on 1/5/26.
//

#ifndef TAIYI_CONTACTS_H
#define TAIYI_CONTACTS_H

#include <cstdint>
#include <vector>
#include <cassert>
#include "Types.h"

//------------------------------------------------------------------------------
// AoS view types (NOT used for storage; only for debug/inspection).
//------------------------------------------------------------------------------
/*struct SoftContact {
    VertexID vertex{};
    ShapeID    shape{};
    Vec3       normal_ws{};   // push-out direction (unit length preferred)
    Vec3       point_ws{};    // contact point on shape (world space; optional but useful)
    float      phi{};         // signed distance after radii/thickness; negative = penetration
    int32_t    tid{-1};       // debug/thread id (optional)
};

struct RigidContact {
    ShapeID shape0{-1}, shape1{-1};

    Vec3  point0_bs{};  // contact point in body0 local frame
    Vec3  point1_bs{};  // contact point in body1 local frame
    Vec3  normal_ws{};  // world-space normal (must be consistent convention)
    float phi{};        // signed distance; negative = penetration

    float thickness0{0.f};
    float thickness1{0.f};

    uint64_t feature_pair_key{0}; // for warm-start / contact matching
    uint32_t feature_key{0};
    int32_t  tid{-1};
};

//------------------------------------------------------------------------------
// SoA containers
// Design goals:
//  - Fixed-capacity arrays (vectors sized to capacity) + count
//  - O(1) reset via count=0 (no deallocation)
//  - Deterministic memory layout per field; good for SIMD/GPU-style kernels
//------------------------------------------------------------------------------

struct SoftContactsSoA {
    // Storage (size == capacity)
    std::vector<VertexID>   vertex;
    std::vector<MeshID>     mesh;       // tell us which mesh this vertex belong to
    std::vector<ShapeID>    shape;
    std::vector<Vec3>       normal_ws;
    std::vector<Vec3>       point_ws;
    std::vector<float>      phi;
    std::vector<int32_t>    tid;

    // Number of valid contacts currently stored (<= capacity)
    int32_t count{0};

    [[nodiscard]] int32_t capacity() const noexcept {
        return static_cast<int32_t>(phi.size());
    }

    void set_capacity(const int32_t cap) {
        assert(cap >= 0);
        vertex.resize(static_cast<size_t>(cap));
        shape.resize(static_cast<size_t>(cap));
        normal_ws.resize(static_cast<size_t>(cap));
        point_ws.resize(static_cast<size_t>(cap));
        phi.resize(static_cast<size_t>(cap));
        tid.resize(static_cast<size_t>(cap));
        count = 0;
    }

    void clear() noexcept { count = 0; }

    [[nodiscard]] bool empty() const noexcept { return count == 0; }

    // Append one contact. Returns false if capacity exceeded.
    bool push(const VertexID p, const ShapeID s, const Vec3& n_ws, const Vec3& pt_ws, const float phi_val, const int32_t tid_val = -1) {
        const int32_t i = count;
        if (i >= capacity()) return false;
        write(i, p, s, n_ws, pt_ws, phi_val, tid_val);
        count = i + 1;
        return true;
    }

    // Write by index (useful if you later do parallel writes with an atomic counter).
    void write(const int32_t i,
               const VertexID p, const ShapeID s,
               const Vec3& n_ws, const Vec3& pt_ws,
               const float phi_val, const int32_t tid_val = -1) {
        assert(i >= 0 && i < capacity());
        vertex[static_cast<size_t>(i)] = p;
        shape[static_cast<size_t>(i)]    = s;
        normal_ws[static_cast<size_t>(i)] = n_ws;
        point_ws[static_cast<size_t>(i)]  = pt_ws;
        phi[static_cast<size_t>(i)]       = phi_val;
        tid[static_cast<size_t>(i)]       = tid_val;
    }

    // Debug: pack the i-th contact into an AoS record.
    [[nodiscard]] SoftContact get(const int32_t i) const {
        assert(i >= 0 && i < count);
        SoftContact c;
        c.vertex  = vertex[static_cast<size_t>(i)];
        c.shape     = shape[static_cast<size_t>(i)];
        c.normal_ws = normal_ws[static_cast<size_t>(i)];
        c.point_ws  = point_ws[static_cast<size_t>(i)];
        c.phi       = phi[static_cast<size_t>(i)];
        c.tid       = tid[static_cast<size_t>(i)];
        return c;
    }
};

struct RigidContactsSoA {
    std::vector<ShapeID> shape0;
    std::vector<ShapeID> shape1;

    std::vector<Vec3>  point0_bs;
    std::vector<Vec3>  point1_bs;
    std::vector<Vec3>  normal_ws;
    std::vector<float> phi;

    std::vector<float>   thickness0;
    std::vector<float>   thickness1;
    std::vector<uint64_t> feature_pair_key;
    std::vector<uint32_t> feature_key;
    std::vector<int32_t>  tid;

    int32_t count{0};

    [[nodiscard]] int32_t capacity() const noexcept {
        return static_cast<int32_t>(phi.size());
    }

    void set_capacity(const int32_t cap) {
        assert(cap >= 0);
        const size_t n = static_cast<size_t>(cap);
        shape0.resize(n);
        shape1.resize(n);
        point0_bs.resize(n);
        point1_bs.resize(n);
        normal_ws.resize(n);
        phi.resize(n);
        thickness0.resize(n);
        thickness1.resize(n);
        feature_pair_key.resize(n);
        feature_key.resize(n);
        tid.resize(n);
        count = 0;
    }

    void clear() noexcept { count = 0; }
    [[nodiscard]] bool empty() const noexcept { return count == 0; }

    bool push(const ShapeID s0, const ShapeID s1,
              const Vec3& p0_bs, const Vec3& p1_bs,
              const Vec3& n_ws, const float phi_val,
              const float th0 = 0.f, const float th1 = 0.f,
              const uint64_t pair_key = 0, const uint32_t key = 0,
              const int32_t tid_val = -1) {
        const int32_t i = count;
        if (i >= capacity()) return false;
        write(i, s0, s1, p0_bs, p1_bs, n_ws, phi_val, th0, th1, pair_key, key, tid_val);
        count = i + 1;
        return true;
    }

    void write(const int32_t i,
               const ShapeID s0, const ShapeID s1,
               const Vec3& p0_bs, const Vec3& p1_bs,
               const Vec3& n_ws, const float phi_val,
               const float th0 = 0.f, const float th1 = 0.f,
               const uint64_t pair_key = 0, const uint32_t key = 0,
               const int32_t tid_val = -1) {
        assert(i >= 0 && i < capacity());
        const size_t k = static_cast<size_t>(i);
        shape0[k] = s0;
        shape1[k] = s1;
        point0_bs[k] = p0_bs;
        point1_bs[k] = p1_bs;
        normal_ws[k] = n_ws;
        phi[k] = phi_val;

        thickness0[k] = th0;
        thickness1[k] = th1;
        feature_pair_key[k] = pair_key;
        feature_key[k] = key;
        tid[k] = tid_val;
    }

    [[nodiscard]] RigidContact get(const int32_t i) const {
        assert(i >= 0 && i < count);
        const size_t k = static_cast<size_t>(i);
        RigidContact c;
        c.shape0 = shape0[k];
        c.shape1 = shape1[k];
        c.point0_bs = point0_bs[k];
        c.point1_bs = point1_bs[k];
        c.normal_ws = normal_ws[k];
        c.phi = phi[k];
        c.thickness0 = thickness0[k];
        c.thickness1 = thickness1[k];
        c.feature_pair_key = feature_pair_key[k];
        c.feature_key = feature_key[k];
        c.tid = tid[k];
        return c;
    }
};

//------------------------------------------------------------------------------
// Unified container
//------------------------------------------------------------------------------
struct Contacts {
    SoftContactsSoA  soft;
    RigidContactsSoA rigid;

    // Keep capacity; reset counts only.
    void Clear() noexcept {
        soft.clear();
        rigid.clear();
    }

    void SetCapacity(const int32_t soft_cap, const int32_t rigid_cap) {
        soft.set_capacity(soft_cap);
        rigid.set_capacity(rigid_cap);
    }
};*/

#endif //TAIYI_CONTACTS_H