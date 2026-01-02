#include <raylib.h>
#include "Application.h"
#include "RenderHelper.hpp"

#include "Sample.h"

inline Vector3 ToRayVec(Vec3& v_pos) {
    return Vector3{v_pos.x(), v_pos.y(), v_pos.z()};
};

int main(){
    constexpr Application::Desc Desc;
    Application app(Desc);
    RegisterAllSamples(app.Registry());
    app.Run(SampleId::DUMMY_SAMPLE);

    return 0;
}
