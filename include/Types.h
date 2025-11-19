//
// Created by 徐天焱 on 2025/11/5.
//

#ifndef TYPES_H
#define TYPES_H

#include <memory>
#include <Eigen/Dense>
#include <vector>

struct mesh_on_cpu;
using IndexBuffer = std::vector<uint16_t>;
using VertexId = uint16_t;
constexpr VertexId INVALID_VERTEX_ID = 0xFFFF;

using Vec3 = Eigen::Vector3f;
using MeshPtr = std::unique_ptr<mesh_on_cpu>;


#endif //TYPES_H
