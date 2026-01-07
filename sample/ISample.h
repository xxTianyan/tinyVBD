//
// Created by xumiz on 2026/1/1.
//

#ifndef TAIYI_ISAMPLE_H
#define TAIYI_ISAMPLE_H

struct AppContext;

class ISample {
public:
    virtual ~ISample() = default;

    // each sample should have a name for ui
    [[nodiscard]] virtual const char* Name() const = 0;

    // call back for starting the sample
    virtual void OnEnter(AppContext& ctx) = 0;

    // call back for ending the sample
    virtual void OnExit(AppContext& ctx) = 0;

    // simulation step
    virtual void Update(AppContext& ctx) = 0;

    // raylib rendering staff
    virtual void Render(AppContext& ctx) = 0;

    // unique ui panel
    virtual void DrawUI(AppContext& ctx) = 0;

    // reload sample
    virtual void Reset(AppContext& ctx) = 0;
};

using SamplePtr = std::unique_ptr<ISample>;

#endif //TAIYI_ISAMPLE_H