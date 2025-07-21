// textureFragmentShader.glsl
#version 150
in vec2 tc;
out vec4 fragColor;

uniform sampler2D tex;

void main()
{
    fragColor = texture(tex, tc);
}
