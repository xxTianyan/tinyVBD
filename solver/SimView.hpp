//
// Created by tianyan on 1/7/26.
//

#ifndef TAIYI_SIMVIEW_HPP
#define TAIYI_SIMVIEW_HPP

#include <span>
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <cassert>

struct Vec3;
struct tetrahedron;
struct triangle;
struct edge;
struct MeshInfo;

struct State;
struct MModel;

namespace sim {

// ---------------------------
// Small utilities
// ---------------------------
inline void throw_if(bool cond, const char* msg) {
    if (cond) throw std::runtime_error(msg);
}

// ---------------------------
// Topology view (read-only)
// ---------------------------
struct TopologyView {
    std::span<const MeshInfo>     mesh_infos;
    std::span<const tetrahedron>  tets;
    std::span<const triangle>     tris;
    std::span<const edge>         edges;

    uint64_t topology_version = 0;
};

// ---------------------------
// Particle views
// ---------------------------

// Read-only particle view
struct ParticleConstView {
    size_t n = 0;

    std::span<const Vec3>  pos;
    std::span<const Vec3>  vel;
    std::span<const Vec3>  force;

    std::span<const float> inv_mass;
};

// Mutable particle view (typical solver input/output)
struct ParticleView {
    size_t n = 0;

    std::span<Vec3>        pos;
    std::span<Vec3>        vel;
    std::span<Vec3>        force;

    std::span<const float> inv_mass;
};

// ---------------------------
// Stage views (pipeline-level)
// ---------------------------

// Forward step usually needs "prev pos" and optionally an "inertia / predicted pos" buffer.
// These buffers are typically solver scratch (Newton's `self.*`).
struct ForwardStepView {
    ParticleView p;

    std::span<Vec3> pos_prev;    // q_prev
    std::span<Vec3> inertia_pos; // predicted / inertia position (optional but common)
};

// You can extend similarly for SolveView / PostStepView when you get there.
// struct SolveView { ParticleView p; ... };

// ---------------------------
// Solver scratch that binds to model topology/particle count
// ---------------------------
struct SolverScratch {
    // Buffers frequently used across stages; put solver-private temporaries here.
    std::vector<Vec3> particle_pos_prev;
    std::vector<Vec3> particle_inertia_pos;

    // Binding metadata: used to know when to resize/rebind.
    uint64_t bound_topology_version = UINT64_MAX;
    size_t   bound_particle_count   = 0;

    // Ensure scratch buffers match current model (particle count + topology version).
    // Call this once per step (or in MakeForwardStepView).
    void BindTo(const MModel& model);
};

// ---------------------------
// MakeView helpers
// ---------------------------

// Topology view from model (read-only)
TopologyView MakeTopologyView(const MModel& model);

// Particle views
ParticleConstView MakeParticleConstView(const MModel& model, const State& state);
ParticleView      MakeParticleView(const MModel& model, State& state);

// Stage views
ForwardStepView   MakeForwardStepView(const MModel& model, State& state, SolverScratch& scratch);

// ============================================================================
// Inline implementations
// ============================================================================

// NOTE: we implement these inline in the header for convenience.
// If you prefer, you can move implementations to a .cpp.

inline TopologyView MakeTopologyView(const MModel& model) {
    TopologyView v;
    v.mesh_infos = model.mesh_infos;
    v.tets       = model.tets;
    v.tris       = model.tris;








#endif //TAIYI_SIMVIEW_HPP