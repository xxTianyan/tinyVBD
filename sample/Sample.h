//
// Created by xumiz on 2025/12/24.
//

#ifndef TINYVBD_SAMPLE_H
#define TINYVBD_SAMPLE_H

#include <raylib.h>
#include "World.h"
#include "ShaderManager.h"
#include "VBDDynamics.h"

class Sample {
public:
    Sample();
    virtual ~Sample() = default;

    virtual void CreateWorld() {};
    virtual void CreateFloor();
    virtual void Step(float dt);
    virtual void CleanUp();

    // for simulation
    std::unique_ptr<World> m_world;
    std::unique_ptr<VBDSolver> m_solver;

    // for rendering
    std::vector<Model> m_models;
    std::unique_ptr<ShaderManager> m_shader_manager;

    bool isPaused = true;
    Model m_floor{};
};


#endif //TINYVBD_SAMPLE_H