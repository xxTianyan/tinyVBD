#version 330

// Raylib mesh attributes
in vec3 vertexPosition;
in vec3 vertexNormal;
in vec2 vertexTexCoord;
in vec4 vertexColor;

// ----------- Shared matrices (Raylib will feed these if shader.locs are set) -----------
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matNormal;

out vec3 vWorldPos;
out vec3 vWorldNormal;
out vec2 vTexCoord;
out vec4 vColor;

void main()
{
    vec4 wp = matModel * vec4(vertexPosition, 1.0);
    vWorldPos = wp.xyz;

    vWorldNormal = normalize((matNormal * vec4(vertexNormal, 0.0)).xyz);

    vTexCoord = vertexTexCoord;
    vColor = vertexColor;

    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
