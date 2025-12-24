#ifndef SHADER_MANAGER_H
#define SHADER_MANAGER_H

#include <raylib.h>

struct ClothShader
{
    Shader shader{};
    int lightDirLoc{-1};
};

struct SceneShaders
{
    ClothShader cloth;
};

SceneShaders LoadSceneShaders();
void UpdateClothLighting(const ClothShader &cloth, const Vector3 &lightDir, const Vector3 &viewPos);
void UnloadSceneShaders(SceneShaders &shaders);

#endif // SHADER_MANAGER_H
