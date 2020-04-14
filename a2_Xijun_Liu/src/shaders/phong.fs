#version 330 core


#define M_PI       3.14159265358979323846f
#define INV_PI     0.31830988618379067154f
#define INV_TWOPI  0.15915494309189533577f

uniform vec3 camPos;
uniform vec3 lightPos;
uniform vec3 lightIntensity;
uniform vec3 rho_d;
uniform vec3 rho_s;
uniform float exponent;

in vec3 vNormal;
in vec3 vPos;
out vec3 color;

void main(){
    vec3 wo = camPos - vPos;
    vec3 wi = lightPos - vPos;
    //vec3 wr = reflect(wi, vNormal);
    vec3 wr = 2*normalize(vNormal) * dot(normalize(wi), normalize(vNormal))-normalize(wi);
    float cos_a = dot(normalize(wo), normalize(wr));
    if (cos_a <0){
        cos_a = 0;
    }
    float cos_theta = dot(normalize(wi), normalize(vNormal)) ;
    vec3 intensity = lightIntensity / pow(distance(lightPos, vPos), 2);
    vec3 BRDF = (rho_d * INV_PI) + (rho_s * (exponent + 2.f) * INV_TWOPI * pow(cos_a, exponent));
    color = intensity * BRDF * cos_theta;
}