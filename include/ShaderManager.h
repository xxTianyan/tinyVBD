#ifndef SHADER_MANAGER_H
#define SHADER_MANAGER_H

#include <raylib.h>
#include <string>
#include <unordered_map>
#include <vector>

#include "raymath.h"

struct CommonShaderParams {
    // 1) 背景/雾颜色：决定暗场舞台氛围；fogColor 建议接近背景
    Vector3 skyAmb   = { 0.10f, 0.12f, 0.18f }; // 上方环境：偏冷更高级
    Vector3 gndAmb   = { 0.05f, 0.05f, 0.055f}; // 下方环境：偏灰
    float ambientStrength = 0.04f;              // 环境亮度：越小越“聚光”对比越强

    // 2) 曝光：整体亮度；值越大越亮（配合 HDR）
    float exposure = 1.6f;

    // 3) 弱方向光补光：避免背光面完全死黑；但不要抢聚光灯戏份
    // lightDir 是光线传播方向 (Light->Scene)
    Vector3 lightDir  = Vector3Normalize(Vector3{ -0.4f, -0.6f, -0.2f });
    Vector3 lightColor  = { 0.25f, 0.25f, 0.25f };   // 很弱！如果你想更暗，可降到 0.1

    // 4) 聚光灯
    Vector3 spotPos   = { 8.0f, 8.0f, 8.0f };      // 灯的位置：越高光圈越大更柔
    Vector3 spotDir   = Vector3Normalize(Vector3{ -1.0f, -1.0f, -1.0f }); // 朝下
    Vector3 spotColor = { 1.0f, 0.95f, 0.85f };    // 略暖色，舞台灯感觉
    float spotIntensity = 22.0f;                   // 光圈亮度强弱（主要旋钮）
    float spotRange     = 30.0f;                   // 光能覆盖多远（太小会断崖式黑）
    float innerDeg = 12.0f;                        // 内锥角：中心最亮区域大小
    float outerDeg = 20.0f;                        // 外锥角：边缘软过渡宽度
    float spotInnerCos = cosf(innerDeg * DEG2RAD);
    float spotOuterCos = cosf(outerDeg * DEG2RAD);
};

class ShaderManager {
public:
    struct ManagedShader {
        Shader shader{};
        int lightDirLoc{-1};
        int viewPosLoc{-1};
    };

    ManagedShader &LoadShaderProgram(const std::string &name, const char *vertexPath, const char *fragmentPath);

    ManagedShader *Get(const std::string &name);
    bool Has(const std::string &name) const;

    void ApplyShaderToModel(const Model &model, const std::string &name);
    void ApplyShaderToModels(const std::vector<Model> &models, const std::string &name);

    void UpdateViewPos(const std::string &name, const Vector3 &viewPos);
    void UnloadAll();

    static void BindMatrices(const Shader &sh);
    static int CheckSetShaderLocation(const Shader &sh, const char* uniform_val);
    static void SetCommonShaderParams(const Shader &sh);

private:
    ManagedShader &EmplaceOrReplace(const std::string &name);

    std::unordered_map<std::string, ManagedShader> m_shaders;
};

#endif // SHADER_MANAGER_H
