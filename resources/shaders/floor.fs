#version 330

in vec3 vWorldPos;
in vec3 vWorldNormal;

out vec4 fragColor;

// ====================== [Shared / Common Lighting Params] ======================
uniform vec3  viewPos;
uniform float iTime;

uniform vec3  skyAmbientColor;
uniform vec3  groundAmbientColor;
uniform float ambientStrength;

uniform vec3  lightDir;          // 弱方向补光（Light -> Scene）
uniform vec3  lightColor;        // 建议较小

uniform vec3  spotPos;
uniform vec3  spotDir;           // (Light -> Scene)
uniform vec3  spotColor;
uniform float spotIntensity;
uniform float spotRange;
uniform float spotInnerCos;
uniform float spotOuterCos;

uniform float exposure;

// ====================== [Floor Look Params] ======================
uniform float tileScale;         // 网格密度：越大格子越小（3~10）
uniform float lineWidth;         // 网格线宽：0.02~0.06
uniform vec3  baseAColor;        // 地板底色 A
uniform vec3  baseBColor;        // 地板底色 B
uniform vec3  lineColor;         // 网格线颜色

uniform float roughness;         // 0..1：地板建议 0.35~0.70
uniform float bumpStrength;      // 0..1：微起伏强度（0.10~0.35）
uniform float fogDensity;        // 雾强度：0.02~0.08
uniform vec3  fogColor;          // 雾颜色：应接近你的背景色

// ------------------ Helpers ------------------
float hash21(vec2 p)
{
    p = fract(p * vec2(123.34, 456.21));
    p += dot(p, p + 45.32);
    return fract(p.x * p.y);
}

float valueNoise(vec2 p)
{
    vec2 i = floor(p);
    vec2 f = fract(p);
    float a = hash21(i + vec2(0,0));
    float b = hash21(i + vec2(1,0));
    float c = hash21(i + vec2(0,1));
    float d = hash21(i + vec2(1,1));
    vec2 u = f*f*(3.0 - 2.0*f);
    return mix(mix(a,b,u.x), mix(c,d,u.x), u.y);
}

float fbm(vec2 p)
{
    float v = 0.0;
    float a = 0.5;
    for (int i = 0; i < 5; i++)
    {
        v += a * valueNoise(p);
        p *= 2.0;
        a *= 0.5;
    }
    return v;
}

// AA grid line
float gridAA(vec2 uv, float width)
{
    vec2 g = fract(uv) - 0.5;
    vec2 ag = abs(g);
    float d = min(ag.x, ag.y);
    float aa = fwidth(d);
    return 1.0 - smoothstep(width - aa, width + aa, d);
}

vec3 fresnelSchlick(float cosTheta, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

void buildTBN(in vec3 N, out vec3 T, out vec3 B)
{
    vec3 up = (abs(N.y) < 0.999) ? vec3(0,1,0) : vec3(1,0,0);
    T = normalize(cross(up, N));
    B = cross(N, T);
}

float spotAttenuation(vec3 worldPos)
{
    vec3 toP = worldPos - spotPos;               // light -> point
    float dist = length(toP);
    vec3 Sd = normalize(spotDir);

    float cosAngle = dot(normalize(toP), Sd);
    float cone = smoothstep(spotOuterCos, spotInnerCos, cosAngle);

    float rangeAtt = 1.0 - clamp(dist / max(spotRange, 1e-6), 0.0, 1.0);
    rangeAtt = rangeAtt * rangeAtt;

    return cone * rangeAtt;
}

void main()
{
    vec3 N0 = normalize(vWorldNormal);
    vec3 V  = normalize(viewPos - vWorldPos);

    // ----- Procedural floor albedo (grid + subtle noise) -----
    vec2 uv = vWorldPos.xz * tileScale;

    float n = fbm(uv * 0.35 + vec2(0.0, iTime * 0.05));
    vec3 albedo = mix(baseAColor, baseBColor, n);

    float line = gridAA(uv, lineWidth);
    albedo = mix(albedo, lineColor, line * 0.55);

    // ----- Micro bump from height field (cheap but effective) -----
    float h0 = fbm(uv * 0.18 + vec2(1.7, 9.2));
    float eps = 0.15;
    float hx = fbm((uv + vec2(eps, 0.0)) * 0.18 + vec2(1.7, 9.2));
    float hz = fbm((uv + vec2(0.0, eps)) * 0.18 + vec2(1.7, 9.2));

    vec3 T, B;
    buildTBN(N0, T, B);
    vec3 N = normalize(N0 - bumpStrength * ((hx - h0) * T + (hz - h0) * B));

    // ----- Ambient (hemi) -----
    float hemi = clamp(N.y * 0.5 + 0.5, 0.0, 1.0);
    vec3 ambientHemi = mix(groundAmbientColor, skyAmbientColor, hemi);
    vec3 ambient = ambientStrength * ambientHemi * albedo;

    // ----- Weak directional fill (very small) -----
    vec3 Ld = normalize(-lightDir);
    float NdotLd = max(dot(N, Ld), 0.0);
    vec3 dirDiffuse = albedo * NdotLd * lightColor;

    vec3 Hd = normalize(Ld + V);
    float shininess = mix(256.0, 16.0, clamp(roughness, 0.0, 1.0));
    float specPowD = pow(max(dot(N, Hd), 0.0), shininess);

    vec3 F0 = vec3(0.04);
    vec3 Fd = fresnelSchlick(max(dot(V, Hd), 0.0), F0);
    vec3 dirSpec = (0.08) * specPowD * Fd * lightColor * NdotLd; // 0.08: keep directional spec tiny

    // ----- Spotlight (main) -----
    vec3 toP = vWorldPos - spotPos;     // light -> point
    float dist = length(toP);
    vec3 Ls = normalize(-toP);          // point -> light
    float spotAtt = spotAttenuation(vWorldPos);

    float NdotLs = max(dot(N, Ls), 0.0);
    vec3 spotDiffuse = albedo * NdotLs * spotColor * (spotIntensity * spotAtt);

    vec3 Hs = normalize(Ls + V);
    float specPowS = pow(max(dot(N, Hs), 0.0), shininess);
    vec3 Fs = fresnelSchlick(max(dot(V, Hs), 0.0), F0);

    // floor spec can be a bit lower if rough, higher if smooth
    float specK = mix(0.18, 0.06, roughness);
    vec3 spotSpec = specK * specPowS * Fs * spotColor * (spotIntensity * spotAtt) * NdotLs;

    // ----- Combine (HDR linear) -----
    vec3 color = ambient + dirDiffuse + dirSpec + spotDiffuse + spotSpec;

    // ----- Fog: helps “infinite” feel + dark stage look -----
    float viewDist = length(viewPos - vWorldPos);
    float fog = 1.0 - exp(-fogDensity * viewDist);
    fog = clamp(fog, 0.0, 1.0);
    color = mix(color, fogColor, fog);

    // ----- Tonemap + gamma -----
    vec3 mapped = vec3(1.0) - exp(-color * max(exposure, 0.0001));
    mapped = pow(mapped, vec3(1.0/2.2));
    fragColor = vec4(mapped, 1.0);
}
