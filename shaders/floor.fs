#version 330

in vec3 vWorldPos;
in vec3 vWorldNormal;

uniform vec3 viewPos;

// 方向光：lightDir = 光线传播方向(光->场景)，所以表面->光源方向 L = -lightDir
uniform vec3  lightDir;
uniform vec3  lightColor;           // 允许 >1 (HDR)

uniform vec3  skyAmbientColor;
uniform vec3  groundAmbientColor;
uniform float ambientStrength;

uniform float iTime;

// 外观控制
uniform float tileScale;            // 例如 3~8
uniform float lineWidth;            // 例如 0.02~0.06（越大线越粗）
uniform float roughness;            // 0~1：越大越“哑”，建议 0.35~0.7
uniform float bumpStrength;         // 0~1：微起伏强度，建议 0.1~0.35
uniform float fogDensity;           // 例如 0.02~0.06
uniform float exposure;             // 1.2~2.2

out vec4 fragColor;

// ---------- util ----------
float hash21(vec2 p)
{
    // cheap hash
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

// 抗锯齿网格线（世界坐标 xz 当 uv）
float gridAA(vec2 uv, float width)
{
    vec2 g = fract(uv) - 0.5;
    vec2 ag = abs(g);
    float d = min(ag.x, ag.y);

    // fwidth 做 AA
    float aa = fwidth(d);
    return 1.0 - smoothstep(width - aa, width + aa, d);
}

vec3 fresnelSchlick(float cosTheta, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

// 在任意法线 N 下构造稳定切线/副切线
void buildTBN(in vec3 N, out vec3 T, out vec3 B)
{
    vec3 up = (abs(N.y) < 0.999) ? vec3(0,1,0) : vec3(1,0,0);
    T = normalize(cross(up, N));
    B = cross(N, T);
}

void main()
{
    vec3 N0 = normalize(vWorldNormal);
    vec3 V  = normalize(viewPos - vWorldPos);
    vec3 L  = normalize(-lightDir);
    vec3 H  = normalize(L + V);

    // ---------- Shadertoy-ish procedural surface ----------
    // world xz -> uv
    vec2 uv = vWorldPos.xz * tileScale;

    // base albedo：暗灰蓝基底 + 轻微噪声变化
    float n = fbm(uv * 0.35 + vec2(0.0, iTime * 0.05));
    vec3 baseA = vec3(0.09, 0.10, 0.11);
    vec3 baseB = vec3(0.13, 0.14, 0.15);
    vec3 albedo = mix(baseA, baseB, n);

    // grid lines（略亮）
    float line = gridAA(uv, lineWidth);
    vec3 lineCol = vec3(0.20, 0.21, 0.23);
    albedo = mix(albedo, lineCol, line * 0.55);

    // ---------- bump (fake micro height -> normal perturb) ----------
    // height field
    float h0 = fbm(uv * 0.18 + vec2(1.7, 9.2));
    // screen-space epsilon (world scale)
    float eps = 0.15;
    float hx = fbm((uv + vec2(eps, 0.0)) * 0.18 + vec2(1.7, 9.2));
    float hz = fbm((uv + vec2(0.0, eps)) * 0.18 + vec2(1.7, 9.2));

    vec3 T, B;
    buildTBN(N0, T, B);

    // gradient -> perturb normal
    float dhx = (hx - h0);
    float dhz = (hz - h0);
    vec3 N = normalize(N0 - bumpStrength * (dhx * T + dhz * B));

    // ---------- lighting ----------
    // hemi ambient
    float hemi = clamp(N.y * 0.5 + 0.5, 0.0, 1.0);
    vec3 ambientHemi = mix(groundAmbientColor, skyAmbientColor, hemi);
    vec3 ambient = ambientStrength * ambientHemi * albedo;

    // diffuse
    float NdotL = max(dot(N, L), 0.0);
    vec3 diffuse = albedo * NdotL * lightColor;

    // specular (roughness -> shininess mapping)
    float shininess = mix(256.0, 16.0, clamp(roughness, 0.0, 1.0)); // rough high -> broad highlight
    float NdotH = max(dot(N, H), 0.0);
    float specPow = pow(NdotH, shininess);

    // dielectric F0
    vec3 F0 = vec3(0.04);
    vec3 F  = fresnelSchlick(max(dot(V, H), 0.0), F0);

    // spec strength decreases with roughness
    float specK = mix(0.25, 0.06, roughness);
    vec3 specular = specK * specPow * F * lightColor * NdotL;

    // ---------- fake reflection of sky (Shadertoy-ish) ----------
    vec3 R = reflect(-V, N);
    float skyT = clamp(R.y * 0.5 + 0.5, 0.0, 1.0);
    vec3 skyCol = mix(groundAmbientColor, skyAmbientColor, smoothstep(0.0, 1.0, skyT));
    // reflection stronger at grazing angles, weaker if rough
    float fres = fresnelSchlick(max(dot(N, V), 0.0), vec3(0.02)).x;
    float reflK = fres * (1.0 - roughness) * 0.65;
    vec3 reflection = skyCol * reflK;

    vec3 color = ambient + diffuse + specular + reflection;

    // ---------- fog (distance fade to background) ----------
    float dist = length(viewPos - vWorldPos);
    float fog = 1.0 - exp(-fogDensity * dist);
    vec3 bg = mix(vec3(0.03,0.04,0.05), vec3(0.10,0.13,0.17), 0.65); // 深色背景基调
    color = mix(color, bg, clamp(fog, 0.0, 1.0));

    // ---------- tonemap + gamma ----------
    vec3 mapped = vec3(1.0) - exp(-color * max(exposure, 0.0001));
    mapped = pow(mapped, vec3(1.0/2.2));
    fragColor = vec4(mapped, 1.0);
}
