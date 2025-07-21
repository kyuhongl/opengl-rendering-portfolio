#version 150

in vec3 position;
in vec2 texCoord;
out vec3 texCoords;

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;

void main()
{
    texCoords = position;
    mat4 viewMatrix = mat4(mat3(modelViewMatrix));
    vec4 pos = projectionMatrix * viewMatrix * vec4(position, 1.0);
    gl_Position = pos.xyww;
} 