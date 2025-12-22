//
// Created by tianyan on 12/22/25.
//

#include "VBDDynamics.h"

/*NodeTetAdj VBDSolver::BuildNodeTetAdj(const size_t num_nodes, const std::vector<tetrahedron>& tets) {

    NodeTetAdj adj;
    adj.offsets.assign(num_nodes+1, 0);

    for (const auto & tet : tets) {
        for (size_t k = 0; k < 4; k++) {
            const VertexId v = tet.vertices[k];
            adj.offsets[v+1]++;
        }
    }

    for (size_t i = 1; i < num_nodes+1; i++) {
        adj.offsets[i] += adj.offsets[i-1];
    }

    adj.incidentTets.resize(adj.offsets.back());

    auto cursor = adj.offsets;

    for (size_t i = 0; i < tets.size(); i++) {
        for (size_t k = 0; k < 4; k++) {
            const VertexId v = tets[i].vertices[k];
            const uint32_t c = cursor[v]++;
            adj.incidentTets[c] = static_cast<uint32_t>(i);
        }
    }
    return adj;
}*/