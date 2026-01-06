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
    static void BuildCloth(float width, float height, int resX, int resY, const Vec3& center = Vec3(0,0,0),
                        ClothOrientation orientation = ClothOrientation::Vertical);

    // build 3D box
    static void BuildBox(float w, float h, float d);

    // Build a UV sphere
    static void BuildSphere(float radius, int sectors, int stacks);

    static void BuildAdjacency(mesh_on_cpu& mesh);

private:
    // static void ApplyFixConsition();

    // static void DistributeMass(mesh_on_cpu& mesh);

};



#endif //TAIYI_MESHBUILDER_H