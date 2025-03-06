#version 150

in vec4 col;
in vec2 texCoordV;
out vec4 c;

uniform float colorCycle; 

// textures
uniform sampler2D textureSampler;
uniform bool useTexture;
uniform bool rainbowMode;

// change hsv to rgb values
vec3 hsv_rgb(vec3 c) {
    vec4 newVec = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    // this is how to change hue values to rgb
    vec3 p = abs(fract(c.xxx + newVec.xyz) * 6.0 - newVec.www);
    return c.z * mix(newVec.xxx, clamp(p - newVec.xxx, 0.0, 1.0), c.y);
}

void main()
{
    if (rainbowMode) { 
        // create rainbow effect by using screen y-coordinate to vary hue, not actual points/pixels
        // colorCycle allows animation of the rainbow pattern over time
        vec3 hsv = vec3(colorCycle + gl_FragCoord.y * 0.001, 1.0, 1.0);
        vec3 rgb = hsv_rgb(hsv);
        c = vec4(rgb, 1.0);
    }
    else if (useTexture) {
        c = texture(textureSampler, texCoordV);  // Using renamed sampler
    }
    else {
        c = col;
    }
}

