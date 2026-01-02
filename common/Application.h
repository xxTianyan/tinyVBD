//
// Created by xumiz on 2026/1/1.
//

#ifndef TINYVBD_APPLICATION_H
#define TINYVBD_APPLICATION_H


#include <memory>
#include "CameraController.h"
#include "SampleRegistry.hpp"
#include "ShaderManager.h"

struct AppContext {
    float dt = 0.0f;
    bool paused = false;

    OrbitCamera orbitCam{};

    ShaderManager* shader_manager = nullptr;
};


class Application {
public:
    struct Desc {
        int width  = 1280;
        int height = 720;
        const char* title = "tinyVBD";
        int target_fps = 120;
        bool resizable = true;
    };

    explicit Application(Desc desc);
    ~Application();

    // get registry sample factory
    SampleRegistry& Registry() { return registry_; }

    // start app
    void Run(SampleId start_sample);

private:
    enum class PendingActionType {
        None,
        SwitchSample,
        ReloadSample,
        ReloadShaders,
        TogglePause,
    };

    struct PendingAction {
        PendingActionType type = PendingActionType::None;
        SampleId target = SampleId::empty_scene; // for switch
    };

private:
    void InitWindowAndUI_() const;
    static void ShutdownWindowAndUI_();

    void EnterSample_(SampleId id);
    void ExitCurrentSample_();

    void RequestSwitchSample_(SampleId id);
    void RequestReloadSample_();

    void PollHotkeys_();
    void DrawAppUI_();     // global ui
    void ExecutePending_(); // delay execute

private:
    Desc desc_;
    SampleRegistry registry_;

    // global resource
    ShaderManager shader_manager_;

    // running state
    AppContext ctx_{};

    SampleId current_id_ = SampleId::empty_scene;
    std::unique_ptr<ISample> current_;

    PendingAction pending_{};
};




#endif //TINYVBD_APPLICATION_H