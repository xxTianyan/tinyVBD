#version 330

// Raylib 默认传入的属性
in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec3 vertexNormal;

// 输出到 FS
out vec3 FragPos;
out vec3 Normal;
out vec2 TexCoord;
out vec4 FragPosLightSpace;

// Uniforms
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 lightVP; // 光源的 View-Projection 矩阵

void main()
{
    // 计算世界坐标
    FragPos = vec3(matModel * vec4(vertexPosition, 1.0));

    // 传递法线
    Normal = normalize(vec3(matModel * vec4(vertexNormal, 0.0)));

    // 传递纹理坐标
    TexCoord = vertexTexCoord;

    // 计算在光源视角下的位置 (用于 FS 里的 ShadowCalculation)
    FragPosLightSpace = lightVP * vec4(FragPos, 1.0);

    // 【只有在这里】设置屏幕裁剪坐标
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}