#ifndef CAMERA_CONTROLLER_H
#define CAMERA_CONTROLLER_H

#include <raylib.h>
#include <raymath.h>
#include <vector>

struct OrbitCtrl {
    float yaw;
    float pitch;
    float radius;
    float rotSens;
    float zoomSens;
    float panSens;
};

struct OrbitCamera {
    Camera3D camera{};
    OrbitCtrl orbit{};
};

OrbitCamera CreateOrbitCamera(Vector3 position, Vector3 target);
void UpdateOrbitCameraMouse(OrbitCamera &orbitCam, Vector2 mouseDelta, float wheelMove,
                            bool rightButtonDown, bool middleButtonDown);
void UpdateOrbitCameraKeyboard(OrbitCamera &orbitCam, float deltaTime, float moveSpeed);
void ReframeOrbitToModels(OrbitCamera &orbitCam, const std::vector<Model> &models, float margin);
void RefreshCameraTransform(OrbitCamera &orbitCam);

#endif // CAMERA_CONTROLLER_H
