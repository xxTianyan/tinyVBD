#version 330

in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec3 vertexNormal;

out vec3 FragPos;
out vec3 Normal;
out vec2 TexCoord;

uniform mat4 mvp;
uniform mat4 matModel;

void main()
{
    FragPos = vec3(matModel * vec4(vertexPosition, 1.0));
    Normal  = normalize(vec3(matModel * vec4(vertexNormal, 0.0)));
    TexCoord = vertexTexCoord;

    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
