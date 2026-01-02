//
// Created by tianyan on 12/23/25.
//

#ifndef MESHBUILDER_H
#define MESHBUILDER_H

#include "Mesh.h"

enum class ClothOrientation {
    Vertical,   //  (XY Plane)
    Horizontal  //  (XZ Plane)
};

class MeshBuilder {
public:
    // 构造一个矩形布料 (2D 网格)
    static void BuildCloth(mesh_on_cpu* mesh, float width, float height, int resX, int resY,
                            const Vec3& center = Vec3(0,0,0), ClothOrientation orientation = ClothOrientation::Vertical);

    // 构造一个实体立方体 (由四面体或六面体分解而成)
    static void BuildBox(mesh_on_cpu* mesh, float w, float h, float d);

    // 构造一个球体 (UV Sphere)
    static void BuildSphere(mesh_on_cpu* mesh, float radius, int sectors, int stacks);

private:
    // 内部辅助函数：重置并初始化网格容器
    static void PrepareMesh(mesh_on_cpu* mesh, size_t num_nodes);
};



#endif //MESHBUILDER_H