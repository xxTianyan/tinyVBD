//
// Created by 徐天焱 on 2025/11/11.
//

#ifndef VBDDYNAMICS_H
#define VBDDYNAMICS_H

#include "Mesh.h"
#include "Types.h"

class VBDSolver {
    int num_iters;

    explicit VBDSolver(const int num_iters) : num_iters(num_iters) {}
    ~VBDSolver() = default;

    void BuildAdjacency(size_t num_nodes, const std::vector<tetrahedron>& tets);



};

struct Node {
    VertexId idx = INVALID_VERTEX_ID;
};





#endif //VBDDYNAMICS_H
