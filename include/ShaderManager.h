#ifndef SHADER_MANAGER_H
#define SHADER_MANAGER_H

#include <raylib.h>

#include <string>
#include <unordered_map>
#include <vector>

class ShaderManager
{
public:
    struct ManagedShader
    {
        Shader shader{};
        int lightDirLoc{-1};
        int viewPosLoc{-1};
    };

    ManagedShader &LoadShaderProgram(const std::string &name, const char *vertexPath, const char *fragmentPath);
    ManagedShader &LoadClothShader(const std::string &name = "cloth");

    ManagedShader *Get(const std::string &name);
    bool Has(const std::string &name) const;

    void ApplyShaderToModel(Model &model, const std::string &name);
    void ApplyShaderToModels(std::vector<Model> &models, const std::string &name);

    void UpdateLighting(const std::string &name, const Vector3 &lightDir, const Vector3 &viewPos);
    void UnloadAll();

private:
    ManagedShader &EmplaceOrReplace(const std::string &name);

    std::unordered_map<std::string, ManagedShader> m_shaders;
};

#endif // SHADER_MANAGER_H
