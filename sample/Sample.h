//
// Created by xumiz on 2025/12/24.
//

#ifndef TINYVBD_SAMPLE_H
#define TINYVBD_SAMPLE_H

#include <raylib.h>
#include "Scene.h"
#include "ShaderManager.h"
#include "VBDDynamics.h"
#include "Application.h"
#include "ISample.h"


class Sample : public ISample {
public:
    Sample() = default;
    ~Sample() override = default;

    // provide sample name
    [[nodiscard]] const char* Name() const override = 0;

    void OnEnter(AppContext& ctx) override;

    void OnExit(AppContext& ctx) override;

    void Update(AppContext& ctx) override;

    void Render(AppContext& ctx) override;

    void DrawUI(AppContext& ctx) override;

    virtual void CreateFloor(AppContext& ctx);

    // api for samples to use
    virtual void CreateWorld([[maybe_unused]]AppContext& ctx) {};

    virtual void Step([[maybe_unused]]float dt) {}

    // clean cpu resource
    virtual void CleanUp();

protected:

    // upload mesh to get raylib model, tip: upload_all_models()
    virtual void BuildRenderResources() {}

    // clean gpu resource
    virtual void DestroyRenderResources();

public:
    // for simulation
    std::unique_ptr<Scene> scene;
    std::unique_ptr<VBDSolver> solver;

    // for rendering
    std::vector<Model> models;
    Model floor{};

private:
    static bool IsModelValid_(const Model& m);

    static void UnloadModelSafe_(Model& m);

};




#endif //TINYVBD_SAMPLE_H