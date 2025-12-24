#include "CameraController.h"

#include <cmath>

namespace
{
    constexpr float kMinRadius = 0.2f;
    constexpr float kMaxRadius = 100.0f;

    void ClampPitch(float *pitch)
    {
        constexpr float lim = DEG2RAD * 89.0f;
        if (*pitch > lim) *pitch = lim;
        if (*pitch < -lim) *pitch = -lim;
    }

    Vector3 SphericalToCartesian(float r, float yaw, float pitch)
    {
        const float cp = cosf(pitch);
        const float sp = sinf(pitch);
        const float cy = cosf(yaw);
        const float sy = sinf(yaw);
        return (Vector3){r * cp * cy, r * sp, r * cp * sy};
    }
}

OrbitCamera CreateOrbitCamera(Vector3 position, Vector3 target)
{
    const Vector3 toCam = Vector3Subtract(position, target);
    const float radius = Vector3Length(toCam);
    const float yaw = atan2f(toCam.z, toCam.x);
    float pitch = radius > 0.0f ? asinf(toCam.y / radius) : 0.0f;
    ClampPitch(&pitch);

    OrbitCamera orbitCam;
    orbitCam.camera = {0};
    orbitCam.camera.position = position;
    orbitCam.camera.target = target;
    orbitCam.camera.up = Vector3{0.0f, 1.0f, 0.0f};
    orbitCam.camera.fovy = 45.0f;
    orbitCam.camera.projection = CAMERA_PERSPECTIVE;

    orbitCam.orbit = {yaw, pitch, radius > 0.001f ? radius : 3.0f, 0.0035f, 0.12f, 0.0025f};
    return orbitCam;
}

void UpdateOrbitCameraInput(OrbitCamera &orbitCam, Vector2 mouseDelta, float wheelMove,
                            bool rightButtonDown, bool middleButtonDown)
{
    if (rightButtonDown)
    {
        orbitCam.orbit.yaw -= mouseDelta.x * orbitCam.orbit.rotSens;
        orbitCam.orbit.pitch -= mouseDelta.y * orbitCam.orbit.rotSens;
        ClampPitch(&orbitCam.orbit.pitch);
    }

    if (fabsf(wheelMove) > 0.0f)
    {
        const float factor = 1.0f - wheelMove * orbitCam.orbit.zoomSens;
        float newR = orbitCam.orbit.radius * factor;
        if (newR < kMinRadius) newR = kMinRadius;
        if (newR > kMaxRadius) newR = kMaxRadius;
        orbitCam.orbit.radius = newR;
    }

    if (middleButtonDown)
    {
        const Vector3 viewDir = Vector3Normalize(SphericalToCartesian(1.0f, orbitCam.orbit.yaw, orbitCam.orbit.pitch));
        const Vector3 right = Vector3Normalize(Vector3CrossProduct(viewDir, orbitCam.camera.up));
        const Vector3 up = Vector3Normalize(Vector3CrossProduct(right, viewDir));
        const Vector3 delta = Vector3Add(Vector3Scale(right, -mouseDelta.x * orbitCam.orbit.panSens * orbitCam.orbit.radius),
                                         Vector3Scale(up, mouseDelta.y * orbitCam.orbit.panSens * orbitCam.orbit.radius));
        orbitCam.camera.target = Vector3Add(orbitCam.camera.target, delta);
    }
}

void RefreshCameraTransform(OrbitCamera &orbitCam)
{
    orbitCam.camera.position = Vector3Add(orbitCam.camera.target,
                                          SphericalToCartesian(orbitCam.orbit.radius, orbitCam.orbit.yaw, orbitCam.orbit.pitch));
    orbitCam.camera.up = {0, 1, 0};
}

void ReframeOrbitToModels(OrbitCamera &orbitCam, const std::vector<Model> &models, float margin)
{
    if (models.empty()) return;

    BoundingBox box0 = GetModelBoundingBox(models[0]);
    Vector3 min = box0.min;
    Vector3 max = box0.max;

    for (size_t i = 1; i < models.size(); ++i)
    {
        BoundingBox bi = GetModelBoundingBox(models[i]);
        min.x = fminf(min.x, bi.min.x);
        min.y = fminf(min.y, bi.min.y);
        min.z = fminf(min.z, bi.min.z);

        max.x = fmaxf(max.x, bi.max.x);
        max.y = fmaxf(max.y, bi.max.y);
        max.z = fmaxf(max.z, bi.max.z);
    }

    const Vector3 center = {(min.x + max.x) * 0.5f, (min.y + max.y) * 0.5f, (min.z + max.z) * 0.5f};
    const Vector3 diag = Vector3Subtract(max, min);
    const float radiusBox = 0.5f * Vector3Length(diag);
    const float fitDist = radiusBox / tanf(DEG2RAD * orbitCam.camera.fovy * 0.5f);
    const float m = (margin > 1.0f ? margin : 1.15f);

    orbitCam.orbit.radius = fitDist * m;
    orbitCam.camera.target = center;
    RefreshCameraTransform(orbitCam);
}
