#version 150

in vec3 position; 
in vec4 color;   // normal

out vec3 viewPosition; 
out vec3 viewNormal;   
uniform mat4 modelViewMatrix;  
uniform mat3 normalMatrix;     
uniform mat4 projectionMatrix; 

void main()
{
  vec4 viewPosition4 = modelViewMatrix * vec4(position, 1.0f);
  viewPosition = viewPosition4.xyz;
  gl_Position = projectionMatrix * viewPosition4;
  vec3 extractedNormal = normalize(color.rgb);
  viewNormal = normalize(normalMatrix * normalize(color.rgb));
}
