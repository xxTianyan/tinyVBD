#include "ShaderManager.h"

SceneShaders LoadSceneShaders()
{
    SceneShaders shaders{};
    shaders.cloth.shader = LoadShader("../shaders/cloth.vs", "../shaders/cloth.fs");
    shaders.cloth.shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(shaders.cloth.shader, "viewPos");
    shaders.cloth.lightDirLoc = GetShaderLocation(shaders.cloth.shader, "lightDir");
    return shaders;
}

void UpdateClothLighting(const ClothShader &cloth, const Vector3 &lightDir, const Vector3 &viewPos)
{
    SetShaderValue(cloth.shader, cloth.shader.locs[SHADER_LOC_VECTOR_VIEW], &viewPos, SHADER_UNIFORM_VEC3);
    if (cloth.lightDirLoc >= 0)
    {
        SetShaderValue(cloth.shader, cloth.lightDirLoc, &lightDir, SHADER_UNIFORM_VEC3);
    }
}

void UnloadSceneShaders(SceneShaders &shaders)
{
    if (shaders.cloth.shader.id > 0)
    {
        UnloadShader(shaders.cloth.shader);
        shaders.cloth.shader = Shader{};
        shaders.cloth.lightDirLoc = -1;
    }
}
