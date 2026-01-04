#include <raylib.h>
#include "Application.h"
#include "RenderHelper.hpp"


inline Vector3 ToRayVec(Vec3& v_pos) {
    return Vector3{v_pos.x(), v_pos.y(), v_pos.z()};
};

int main(){
    constexpr Application::Desc Desc;
    Application app(Desc);
    RegisterAllSamples(app.Registry());
    app.Run(SampleId::BASIC_CLOTH_EXAMPLE);

    return 0;
}
