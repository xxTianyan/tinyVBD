#version 330

in vec3 vWorldPos;
in vec3 vWorldNormal;
in vec2 vTexCoord;
in vec4 vColor;

uniform vec3 viewPos;                  // camera position in world space

// Directional light
// Definition: lightDir is the light RAY direction (from light -> scene), normalized.
// So the direction from surface -> light is L = -lightDir.
uniform vec3  lightDir;                // normalized
uniform vec3  lightColor;              // linear HDR intensity allowed (can be > 1.0)

// Hemispheric ambient (cheap "environment")
uniform vec3  skyAmbientColor;         // linear
uniform vec3  groundAmbientColor;      // linear
uniform float ambientStrength;         // 0..1, usually 0.03..0.12 for dark scenes

// Material controls
uniform float shininess;               // 32..128 typical
uniform float specStrength;            // 0.15..0.6 for cloth-ish look
uniform float wrapDiffuse;             // 0..0.5, 0.15..0.35 recommended for cloth

// Post / output
uniform float exposure;                // 1.0..2.5 typical

uniform sampler2D texture0;
uniform vec4 colDiffuse;               // raylib supplies from material.maps[ALBEDO].color

out vec4 fragColor;

// Fresnel (Schlick)
vec3 fresnelSchlick(float cosTheta, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

void main()
{
    // Double-sided normal handling
    vec3 N = normalize(gl_FrontFacing ? vWorldNormal : -vWorldNormal);

    // Base color
    vec4 texel = texture(texture0, vTexCoord);
    vec4 base  = colDiffuse * texel;         // includes alpha
    vec3 albedo = base.rgb;

    // View and light
    vec3 V = normalize(viewPos - vWorldPos);
    vec3 L = normalize(-lightDir);           // surface -> light
    vec3 H = normalize(L + V);

    // ---- Diffuse (Lambert + optional wrap) ----
    float NdotL = dot(N, L);
    float diff;
    if (wrapDiffuse > 0.0)
    {
        // "wrap" diffuse for softer cloth shading
        diff = clamp((NdotL + wrapDiffuse) / (1.0 + wrapDiffuse), 0.0, 1.0);
    }
    else
    {
        diff = max(NdotL, 0.0);
    }
    vec3 diffuse = albedo * diff * lightColor;

    // ---- Specular (Blinn-Phong + Fresnel) ----
    float NdotH = max(dot(N, H), 0.0);
    float specPow = pow(NdotH, max(shininess, 1.0));

    // Dielectric F0 ~ 0.02..0.08; cloth usually closer to ~0.04
    vec3 F0 = vec3(0.04);
    float VdotH = max(dot(V, H), 0.0);
    vec3 F = fresnelSchlick(VdotH, F0);

    // Simple energy-ish normalization so shininess changes don't blow intensity too much
    float normFactor = (shininess + 8.0) / 8.0;

    // Gate spec by NdotL to avoid highlights on the dark side
    vec3 specular = specStrength * normFactor * specPow * F * lightColor * max(NdotL, 0.0);

    // ---- Hemispheric Ambient ----
    // Use normal.y to blend sky/ground environment
    float hemi = clamp(N.y * 0.5 + 0.5, 0.0, 1.0);
    vec3 ambientHemi = mix(groundAmbientColor, skyAmbientColor, hemi);
    vec3 ambient = ambientStrength * ambientHemi * albedo;

    // ---- Combine (linear HDR) ----
    vec3 color = ambient + diffuse + specular;

    // ---- Tonemap + Gamma ----
    // Exposure + exponential tonemap (stable for HDR)
    vec3 mapped = vec3(1.0) - exp(-color * max(exposure, 0.0001));
    // Gamma encode to sRGB-ish
    mapped = pow(mapped, vec3(1.0 / 2.2));

    fragColor = vec4(mapped, base.a);
}
