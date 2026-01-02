//
// Created by xumiz on 2025/12/24.
//

#ifndef SAMPLE_H
#define SAMPLE_H

#include <memory>
#include <raylib.h>
#include <vector>
#include "ISample.h"

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
    // TODO: set a good way to init, maybe in detailed exmaples?
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




#endif //SAMPLE_H