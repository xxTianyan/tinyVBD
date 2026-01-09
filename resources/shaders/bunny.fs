#version 330

in vec3 fragPosWS;
in vec3 fragNormalWS;
in vec2 fragUV;
in vec4 fragColor;

out vec4 finalColor;

// ====================== [Raylib Material Inputs] ======================
uniform sampler2D texture0;
uniform vec4 colDiffuse;

// ====================== [Shared / Common Lighting Params] ======================
uniform vec3  viewPos;
uniform float iTime;

uniform vec3  skyAmbientColor;
uniform vec3  groundAmbientColor;
uniform float ambientStrength;

uniform vec3  lightDir;          // Light -> Scene
uniform vec3  lightColor;

uniform vec3  spotPos;
uniform vec3  spotDir;           // Light -> Scene
uniform vec3  spotColor;
uniform float spotIntensity;
uniform float spotRange;
uniform float spotInnerCos;
uniform float spotOuterCos;

uniform float exposure;

// =========================================================
// [TUNING] 你只需要改这一行：
//   0.0 = 更柔和（面片弱）
//   1.0 = 更强面片感（面片明显）
// =========================================================
const float PLASTER_FACET = 0.75;   // 建议：0.25~0.95 之间调

// ---- helpers ----
float saturate(float x) { return clamp(x, 0.0, 1.0); }

float hash13(vec3 p)
{
    p = fract(p * 0.1031);
    p += dot(p, p.yzx + 33.33);
    return fract((p.x + p.y) * p.z);
}

vec3 tonemapExp(vec3 x, float expv)
{
    return vec3(1.0) - exp(-x * max(expv, 1e-4));
}

void main()
{
    vec4 tex = texture(texture0, fragUV);
    vec3 albedo = tex.rgb * colDiffuse.rgb * fragColor.rgb;
    float alpha = tex.a * colDiffuse.a * fragColor.a;

    // 石膏/雕塑更接近“非纯白”，避免向光侧炸白（重要）
    albedo = min(albedo, vec3(0.90));

    vec3 V = normalize(viewPos - fragPosWS);

    // ---------------------------
    // 1) 法线：smooth + 几何面片法线 混合
    // ---------------------------
    vec3 Nv = normalize(fragNormalWS);

    // per-triangle geometric normal -> 面片感来源
    vec3 Ng = normalize(cross(dFdx(fragPosWS), dFdy(fragPosWS)));
    if (!gl_FrontFacing) Ng = -Ng;

    // diffuse 用的法线（可控面片强度）
    float facetStrength = mix(0.18, 0.70, PLASTER_FACET);
    vec3 N = normalize(mix(Nv, Ng, facetStrength));

    // spec 用更“面片”的法线（高光更容易体现面片，但不油）
    float facetSpec = saturate(facetStrength + 0.20);
    vec3 Nspec = normalize(mix(Nv, Ng, facetSpec));

    float NdotV = saturate(dot(N, V));

    // ---------------------------
    // 2) 半球环境光（保证背光不死黑）
    // ---------------------------
    float hemiT = saturate(N.y * 0.5 + 0.5);
    vec3 hemi = mix(groundAmbientColor, skyAmbientColor, hemiT);

    vec3 ambient = hemi * ambientStrength + 0.03 * hemi;

    // ---------------------------
    // 3) 方向光（wrap + 对比 shaping）
    // ---------------------------
    vec3 Ld = normalize(-lightDir);     // scene -> light incoming is -lightDir
    float ndl_d = dot(N, Ld);

    // 越面片（PLASTER_FACET 越大），wrap 越小，面片感更明显
    float wrap = mix(0.30, 0.18, PLASTER_FACET);
    float diff_d = saturate((ndl_d + wrap) / (1.0 + wrap));

    // 过度柔化会糊面：这里收紧一些，并随面片程度调整
    float edge0 = mix(0.02, 0.08, PLASTER_FACET);
    diff_d = smoothstep(edge0, 0.98, diff_d);

    // 对比 shaping（面片强时稍微更“硬”一点）
    float diffPow = mix(1.35, 1.60, PLASTER_FACET);
    diff_d = pow(diff_d, diffPow);

    // 背光抬升，避免纯黑（但别太多）
    float back_d = pow(saturate(-ndl_d), 1.2) * mix(0.22, 0.14, PLASTER_FACET);

    vec3 dirLight = (diff_d + back_d) * lightColor;

    // ---------------------------
    // 4) 聚光灯（cone + range）
    // ---------------------------
    vec3 toLight = spotPos - fragPosWS;
    float dist = length(toLight);
    vec3 Ls = (dist > 1e-6) ? (toLight / dist) : vec3(0.0, 1.0, 0.0);

    float att = saturate(1.0 - dist / max(spotRange, 1e-4));
    att *= att;

    vec3 lightToFrag = normalize(fragPosWS - spotPos);
    float coneCos = dot(lightToFrag, normalize(spotDir));
    float cone = smoothstep(spotOuterCos, spotInnerCos, coneCos);

    float ndl_s = dot(N, Ls);
    float diff_s = saturate((ndl_s + wrap) / (1.0 + wrap));
    diff_s = smoothstep(edge0, 0.98, diff_s);
    diff_s = pow(diff_s, mix(1.30, 1.55, PLASTER_FACET));

    vec3 spotLight = spotColor * (spotIntensity * att * cone * diff_s);

    // ---------------------------
    // 5) 石膏高光：显著减弱 + soft-knee 防炸白
    // ---------------------------
    float fres = pow(1.0 - saturate(dot(normalize(mix(N, Nspec, 0.5)), V)), 5.0);

    vec3 Hd = normalize(Ld + V);
    vec3 Hs = normalize(Ls + V);

    // 关键：把强度拉低到石膏量级；shininess 提高，避免整片发亮
    float specStrength = mix(0.010, 0.020, PLASTER_FACET);
    float shinD = mix(40.0, 60.0, PLASTER_FACET);
    float shinS = mix(44.0, 70.0, PLASTER_FACET);

    float spec_d = pow(saturate(dot(Nspec, Hd)), shinD) * specStrength;
    float spec_s = pow(saturate(dot(Nspec, Hs)), shinS) * (specStrength * 1.2) * (att * cone);

    float spec = (spec_d + spec_s) * (0.40 + 0.70 * fres);

    // soft-knee 压缩：强高光不会直接顶白
    spec = spec / (1.0 + spec);

    // 轮廓光也收敛，防止“发白边”
    float rim = pow(1.0 - NdotV, mix(2.8, 2.2, 1.0 - PLASTER_FACET)) * 0.05;

    // ---------------------------
    // 6) 合成 + 限制黑位（不死黑）+ 细微颗粒
    // ---------------------------
    vec3 lighting = ambient + dirLight + spotLight;

    // 防止压黑过度（尤其背光侧）
    lighting = max(lighting, vec3(0.06));

    vec3 color = albedo * lighting;

    // spec/rim 作为中性（灰白）叠加
    color += vec3(spec);
    color += rim * (0.35 * hemi);

    // plaster grain
    float n = hash13(fragPosWS * 2.0 + vec3(iTime * 0.02));
    color *= 1.0 + (n - 0.5) * 0.03;

    // 强去饱和 -> 黑白灰雕塑感
    float lum = dot(color, vec3(0.299, 0.587, 0.114));
    color = mix(color, vec3(lum), 0.88);

    // exposure + tonemap + gamma
    color = tonemapExp(color, exposure);
    color = pow(color, vec3(1.0 / 2.2));

    finalColor = vec4(color, alpha);
}
