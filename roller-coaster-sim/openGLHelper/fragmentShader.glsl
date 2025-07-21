#version 150

in vec3 viewPosition; // vertex program outputs
in vec3 viewNormal;   // n: vertex program outputs 

out vec4 c; // output color

uniform vec4 La; // light ambient
uniform vec4 Ld; // light diffuse
uniform vec4 Ls; // light specular
uniform vec3 viewLightDirection; 

uniform vec4 ka;
uniform vec4 kd;
uniform vec4 ks;
uniform float alpha;

void main()
{
    // v: camera is at (0,0,0) after the modelview transformation
    vec3 eyedir = normalize(vec3(0,0,0) - viewPosition);

    // r: reflected light direction
    vec3 reflectDir = -reflect(normalize(viewLightDirection), viewNormal);

    // l * n
    float d = max(dot(normalize(viewLightDirection), viewNormal), 0.0f);
    // r * v
    float s = max(dot(reflectDir, eyedir), 0.0f);

    // I: compute the final color
    c = ka * La + d * kd * Ld + pow(s, alpha) * ks * Ls;
}
