//
// Created by xumiz on 2026/1/1.
//

#ifndef TAIYI_APPLICATION_H
#define TAIYI_APPLICATION_H


#include <memory>
#include "CameraController.h"
#include "PefMonitor.hpp"
#include "RenderHelper.h"
#include "SampleRegistry.h"
#include "ShaderManager.h"

struct AppContext {
    float dt = 0.0f;
    bool paused = true;
    int target_fps = 0;

    OrbitCamera* orbitCam = nullptr;

    ShaderManager* shader_manager = nullptr;
};


class Application {
public:
    struct Desc {
        int width  = 1280;
        int height = 720;
        const char* title = "ProjectTaiyi";
        int target_fps = 120;
        bool resizable = true;
    };

    explicit Application(const Desc &desc);
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
        SampleId target = SampleId::EMPTY_SCENE; // for switch
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
    PerformanceMonitor perfMonitor_;
    OrbitCamera orbitCam_;

    SampleId current_id_ = SampleId::EMPTY_SCENE;
    std::unique_ptr<ISample> current_;

    PendingAction pending_{};
};




#endif //TAIYI_APPLICATION_H
