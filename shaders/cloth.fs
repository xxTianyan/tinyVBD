#version 330

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoord;

out vec4 finalColor;

uniform vec3 lightDir;
uniform vec3 viewPos;
uniform vec4 colDiffuse;

void main()
{
    vec3 normal = normalize(gl_FrontFacing ? Normal : -Normal);
    vec3 viewDirection = normalize(viewPos - FragPos);
    vec3 lightDirection = normalize(lightDir);

    vec3 baseColor = colDiffuse.rgb;

    float diffuse = max(dot(normal, lightDirection), 0.0);
    float wrappedDiffuse = clamp(diffuse * 0.75 + 0.25, 0.0, 1.0);

    float fresnel = pow(1.0 - max(dot(normal, viewDirection), 0.0), 3.0);
    float backLight = pow(clamp(dot(-lightDirection, normal), 0.0, 1.0), 2.0);

    vec3 ambient = baseColor * 0.35;
    vec3 clothDiffuse = baseColor * wrappedDiffuse;
    vec3 clothSheen = baseColor * (0.15 * fresnel + 0.1 * backLight);

    vec3 color = ambient + clothDiffuse + clothSheen;
    finalColor = vec4(color, 1.0);
}
