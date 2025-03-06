#version 150

in vec3 position;
in vec4 color;
in vec3 p_center;
in vec3 p_left;
in vec3 p_right;
in vec3 p_up;
in vec3 p_down;

out vec4 col;

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform int mode;
uniform float scale;
uniform float exponent;

// textures
in vec2 texCoord;  
out vec2 texCoordV;

void main()
{
    vec3 pos;
    
    if (mode == 1) {
        // Smoothing mode
        pos = (p_center + p_left + p_right + p_up + p_down) / 5.0;
        
        float height = pos.y;
        height = scale * pow(height, exponent);
        pos.y = height;
        
        // Color based on height
        float brightness = pow(height, exponent);
        col = vec4(brightness, brightness, brightness, 1.0);
    } else {
        pos = position;
        col = color;
    }
    
    gl_Position = projectionMatrix * modelViewMatrix * vec4(pos, 1.0);

    texCoordV = texCoord;
}