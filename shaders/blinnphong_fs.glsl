#version 330 core

in vec3 vFragPos;
in vec3 vNormal;

out vec4 outColor;

// 光照与相机
uniform vec3 uLightPos;     // 点光源（世界）
uniform vec3 uLightColor;   // 0..1
uniform vec3 uViewPos;      // 相机（世界）

// 材质
uniform vec3  uAlbedoColor;   // 纯色漫反射，比如 vec3(1.0,0.7,0.2)
uniform vec3  uAmbient;       // 环境光强度，比如 vec3(0.15)
uniform float uShininess;     // 高光指数，如 64
uniform float uSpecStrength;  // 高光强度，如 0.7
uniform vec4 colDiffuse;  // raylib 从 DrawModel 的 color 注入

void main() {
    vec3 N = normalize(vNormal);
    vec3 L = normalize(uLightPos - vFragPos);
    vec3 V = normalize(uViewPos  - vFragPos);
    vec3 H = normalize(L + V);

    float diff = max(dot(N, L), 0.0);

    float spec = 0.0;
    if (diff > 0.0) {
        spec = pow(max(dot(N, H), 0.0), uShininess) * uSpecStrength;
    }

    float d = length(uLightPos - vFragPos);
    float attenuation = 1.0 / (1.0 + 0.09 * d + 0.032 * d * d);

    // ★ 用 DrawModel 的 tint 调制材质色
    vec3 baseColor = uAlbedoColor * colDiffuse.rgb;

    vec3 ambient  = uAmbient * baseColor;
    vec3 diffuse  = diff * baseColor * uLightColor;
    vec3 specular = spec * uLightColor;

    vec3 color = ambient + (diffuse + specular) * attenuation;
    outColor = vec4(color, 1.0) * colDiffuse.a;  // 可选：用 tint 的 alpha
}

