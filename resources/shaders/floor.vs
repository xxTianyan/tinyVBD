#version 330

in vec3 vertexPosition;
in vec3 vertexNormal;

uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matNormal;

out vec3 vWorldPos;
out vec3 vWorldNormal;

void main()
{
    vec4 wp = matModel * vec4(vertexPosition, 1.0);
    vWorldPos = wp.xyz;
    vWorldNormal = normalize((matNormal * vec4(vertexNormal, 0.0)).xyz);
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
