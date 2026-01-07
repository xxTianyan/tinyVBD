//
// Created by xumiz on 2025/12/24.
//

#ifndef TAIYI_SAMPLE_H
#define TAIYI_SAMPLE_H

#include <memory>
#include <raylib.h>
#include <vector>
#include "ISample.h"
#include "RenderHelper.h"

class Scene;
class VBDSolver;
struct AppContext;

class Sample : public ISample {
public:
    Sample() = default;
    ~Sample() override = default;

    // provide sample name
    [[nodiscard]] const char* Name() const override { return "Empty Scene";};

    void OnEnter(AppContext& ctx) override;

    void OnExit(AppContext& ctx) override;

    void Update(AppContext& ctx) override;

    void Render(AppContext& ctx) override;

    void DrawUI(AppContext& ctx) override;

    void Reset(AppContext& ctx) override;

    virtual void CreateFloor(AppContext& ctx);

    // api functions that need to be over-ride when inherited
    virtual void CreateWorld([[maybe_unused]]AppContext& ctx) {};
    virtual void Step([[maybe_unused]]const float dt) {}
    virtual void BindShaders([[maybe_unused]]AppContext& ctx) {};

    // clean cpu resource
    virtual void CleanUp();

protected:
    virtual void BuildRenderResources();

    // clean gpu resource
    virtual void DestroyRenderResources();


public:
    // for simulation， remember to initialize
    std::unique_ptr<Scene> scene_;
    std::unique_ptr<VBDSolver> solver_;

    // for rendering
    RenderHelper renderHelper_;
    Model floor_{};

};




#endif //TAIYI_SAMPLE_H