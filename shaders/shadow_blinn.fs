#version 330

// 输入：必须与 Vertex Shader 的输出变量名完全一致
in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoord;
in vec4 FragPosLightSpace; // 这是我们在 VS 里算好传过来的，不要用 gl_Position

out vec4 finalColor;

// Uniforms
uniform sampler2D shadowMap;
uniform vec3 lightDir;
uniform vec3 viewPos;
uniform vec4 colDiffuse;

// 阴影计算函数
float ShadowCalculation(vec4 fragPosLightSpace)
{
    // 1. 透视除法 (Perspective Divide)
    vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;

    // 2. 变换到 [0,1] 范围
    projCoords = projCoords * 0.5 + 0.5;

    // 3. 如果超过远平面的深度，认为不在阴影中
    if(projCoords.z > 1.0)
        return 0.0;

    // 4. 计算 Bias (防止阴影波纹/Shadow Acne)
    vec3 normal = normalize(Normal);
    vec3 lightDirection = normalize(lightDir);
    float bias = max(0.005 * (1.0 - dot(normal, lightDirection)), 0.0005);

    // 5. PCF 柔化阴影 (简单的 3x3 采样)
    float shadow = 0.0;
    vec2 texelSize = 1.0 / textureSize(shadowMap, 0);
    float currentDepth = projCoords.z;

    for(int x = -1; x <= 1; ++x)
    {
        for(int y = -1; y <= 1; ++y)
        {
            float pcfDepth = texture(shadowMap, projCoords.xy + vec2(x, y) * texelSize).r;
            shadow += currentDepth - bias > pcfDepth ? 1.0 : 0.0;
        }
    }
    shadow /= 9.0;

    return shadow;
}

void main()
{
    vec3 color = colDiffuse.rgb;
    vec3 normal = normalize(gl_FrontFacing ? Normal : -Normal);
    vec3 lightDirection = normalize(lightDir);

    // 1. Ambient (环境光)
    vec3 ambient = 0.3 * color;

    // 2. Diffuse (漫反射)
    float diff = max(dot(normal, lightDirection), 0.0);
    vec3 diffuse = diff * color;

    // 3. Specular (高光 - Blinn-Phong)
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 halfwayDir = normalize(lightDirection + viewDir);
    float spec = pow(max(dot(normal, halfwayDir), 0.0), 32.0);
    vec3 specular = vec3(0.3) * spec; // 假设高光强度 0.3

    // 4. Shadow (计算阴影因子)
    // 【关键点】这里传入 FragPosLightSpace，绝对不要写成 gl_Position
    float shadow = ShadowCalculation(FragPosLightSpace);

    // 5. 组合结果
    vec3 lighting = (ambient + (1.0 - shadow) * (diffuse + specular));

    finalColor = vec4(lighting, 1.0);
}