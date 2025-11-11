//
// Created by 徐天焱 on 2025/11/5.
//

#ifndef TYPES_H
#define TYPES_H

#include <Eigen/Dense>
#include <vector>

using IndexBuffer = std::vector<uint16_t>;
using VertexId = uint16_t;
constexpr VertexId INVALID_VERTEX_ID = 0xFFFF;

using Vec3 = Eigen::Vector3f;



#endif //TYPES_H
