#version 330 core

layout (location = 0) in vec3 vertexPosition;
layout (location = 2) in vec3 vertexNormal;

uniform mat4 mvp;       // raylib 自动设置
uniform mat4 matModel;  // raylib 自动设置

out vec3 vFragPos;      // 世界坐标
out vec3 vNormal;       // 世界空间法线

void main() {
    vec4 worldPos = matModel * vec4(vertexPosition, 1.0);
    vFragPos = worldPos.xyz;

    // 法线矩阵处理非均匀缩放
    mat3 normalMatrix = transpose(inverse(mat3(matModel)));
    vNormal = normalize(normalMatrix * vertexNormal);

    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
