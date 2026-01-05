//
// Created by tianyan on 12/23/25.
//

#ifndef TAIYI_MESHBUILDER_H
#define TAIYI_MESHBUILDER_H

#include "Types.h"
struct mesh_on_cpu;

enum class ClothOrientation {
    Vertical,   //  (XY Plane)
    Horizontal  //  (XZ Plane)
};

class MeshBuilder {
public:
    // build 2D cloth
    static void BuildCloth(mesh_on_cpu* mesh, float width, float height, int resX, int resY,
                            const Vec3& center = Vec3(0,0,0), ClothOrientation orientation = ClothOrientation::Vertical);

    // build 3D box
    static void BuildBox(mesh_on_cpu* mesh, float w, float h, float d);

    // Build a UV sphere
    static void BuildSphere(mesh_on_cpu* mesh, float radius, int sectors, int stacks);

private:
    // help prepare mesh
    static void PrepareMesh(mesh_on_cpu* mesh, size_t num_nodes);
};



#endif //TAIYI_MESHBUILDER_H