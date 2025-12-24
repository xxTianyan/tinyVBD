#include "ShaderManager.h"

ShaderManager::ManagedShader &ShaderManager::LoadShaderProgram(const std::string &name, const char *vertexPath, const char *fragmentPath)
{
    ManagedShader &managed = EmplaceOrReplace(name);
    managed.shader = LoadShader(vertexPath, fragmentPath);
    managed.viewPosLoc = GetShaderLocation(managed.shader, "viewPos");
    managed.shader.locs[SHADER_LOC_VECTOR_VIEW] = managed.viewPosLoc;
    managed.lightDirLoc = GetShaderLocation(managed.shader, "lightDir");
    return managed;
}

ShaderManager::ManagedShader &ShaderManager::LoadClothShader(const std::string &name)
{
    return LoadShaderProgram(name, "../shaders/cloth.vs", "../shaders/cloth.fs");
}

ShaderManager::ManagedShader *ShaderManager::Get(const std::string &name)
{
    auto it = m_shaders.find(name);
    if (it == m_shaders.end())
    {
        return nullptr;
    }
    return &it->second;
}

bool ShaderManager::Has(const std::string &name) const
{
    return m_shaders.find(name) != m_shaders.end();
}

void ShaderManager::ApplyShaderToModel(Model &model, const std::string &name)
{
    ManagedShader *shader = Get(name);
    if (shader)
    {
        model.materials[0].shader = shader->shader;
    }
}

void ShaderManager::ApplyShaderToModels(std::vector<Model> &models, const std::string &name)
{
    for (auto &model : models)
    {
        ApplyShaderToModel(model, name);
    }
}

void ShaderManager::UpdateLighting(const std::string &name, const Vector3 &lightDir, const Vector3 &viewPos)
{
    ManagedShader *shader = Get(name);
    if (!shader)
    {
        return;
    }

    if (shader->viewPosLoc >= 0)
    {
        SetShaderValue(shader->shader, shader->viewPosLoc, &viewPos, SHADER_UNIFORM_VEC3);
    }

    if (shader->lightDirLoc >= 0)
    {
        SetShaderValue(shader->shader, shader->lightDirLoc, &lightDir, SHADER_UNIFORM_VEC3);
    }
}

void ShaderManager::UnloadAll()
{
    for (auto &entry : m_shaders)
    {
        if (entry.second.shader.id > 0)
        {
            UnloadShader(entry.second.shader);
        }
    }
    m_shaders.clear();
}

ShaderManager::ManagedShader &ShaderManager::EmplaceOrReplace(const std::string &name)
{
    auto [it, inserted] = m_shaders.emplace(name, ManagedShader{});
    if (!inserted)
    {
        if (it->second.shader.id > 0)
        {
            UnloadShader(it->second.shader);
        }
        it->second = ManagedShader{};
    }
    return it->second;
}
