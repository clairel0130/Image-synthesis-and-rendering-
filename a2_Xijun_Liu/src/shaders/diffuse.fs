#version 330 core


#define M_PI       3.14159265358979323846f
#define INV_PI     0.31830988618379067154f
#define INV_TWOPI  0.15915494309189533577f

uniform vec3 camPos;
uniform vec3 lightPos;
uniform vec3 lightIntensity;
uniform vec3 albedo;

in vec3 vNormal;
in vec3 vPos;
out vec3 color;

void main(){
    float cos_theta = dot((lightPos - vPos), vNormal) / (length(lightPos - vPos) * length(vNormal));
    vec3 intensity = lightIntensity / (distance(lightPos, vPos)*distance(lightPos, vPos));
    color = cos_theta * albedo * intensity / M_PI;
}