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
using VertexID = uint16_t;
using MeshID = uint32_t;
using MaterialID = uint32_t;
constexpr VertexID INVALID_VERTEX_ID = 0xFFFF;
constexpr uint32_t INVALID_INDEX = 0xFFFFFFFF;

using Vec3 = Eigen::Vector3f;
using Mat3 = Eigen::Matrix3f;
using Mat2 = Eigen::Matrix2f;
using Mat32 = Eigen::Matrix<float, 3, 2>;
using MeshPtr = std::unique_ptr<mesh_on_cpu>;



#endif //TYPES_H
