#version 150

in vec3 texCoords;
out vec4 fragColor;

uniform samplerCube skyboxTexture;

void main()
{
    vec3 adjustedCoords = texCoords;
    
    if (abs(texCoords.y) > abs(texCoords.x) && abs(texCoords.y) > abs(texCoords.z) && texCoords.y > 0) {
        adjustedCoords.x = -texCoords.z;
        adjustedCoords.z = texCoords.x;
    }
    else if (abs(texCoords.y) > abs(texCoords.x) && abs(texCoords.y) > abs(texCoords.z) && texCoords.y < 0) {
        adjustedCoords.x = -texCoords.z;
        adjustedCoords.z = texCoords.x;
    }
    
    vec4 texColor = texture(skyboxTexture, adjustedCoords);
    
    if(texColor.r == 0.0 && texColor.g == 0.0 && texColor.b == 0.0) {
        fragColor = vec4(1.0, 0.0, 0.0, 1.0); // bright red for debug
    } else {
        fragColor = texColor;
    }
} 