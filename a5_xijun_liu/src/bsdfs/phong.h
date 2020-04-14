/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Modified Phong reflectance model
 */
struct PhongBSDF : BSDF {

    std::unique_ptr<Texture < v3f>> specularReflectance;
    std::unique_ptr<Texture < v3f>> diffuseReflectance;
    std::unique_ptr<Texture < float>> exponent;
    float specularSamplingWeight;
    float scale;

    PhongBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
        const tinyobj::material_t& mat = scene.materials[matID];

        if (mat.specular_texname.empty())
            specularReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.specular)));
        else
            specularReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.specular_texname));

        if (mat.diffuse_texname.empty())
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.diffuse)));
        else
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.diffuse_texname));

        exponent = std::unique_ptr<Texture<float>>(new ConstantTexture1f(mat.shininess));

        //get scale value to ensure energy conservation
        v3f maxValue = specularReflectance->getMax() + diffuseReflectance->getMax();
        float actualMax = max(max(maxValue.x, maxValue.y), maxValue.z);
        scale = actualMax > 1.0f ? 0.99f * (1.0f / actualMax) : 1.0f;

        float dAvg = getLuminance(diffuseReflectance->getAverage() * scale);
        float sAvg = getLuminance(specularReflectance->getAverage() * scale);
        specularSamplingWeight = sAvg / (dAvg + sAvg);

        components.push_back(EGlossyReflection);
        components.push_back(EDiffuseReflection);

        combinedType = 0;
        for (unsigned int component : components)
            combinedType |= component;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    v3f eval(const SurfaceInteraction& i) const override {
        v3f val(0.f);
        // TODO: Add previous assignment code (if needed)
        float exp = exponent->eval(worldData, i);
        float cosA = 0.f;
        float cos = 0.f;
        // TODO: Implement this
        v3f wo(0.f);
        if ((Frame::cosTheta(i.wi)>= 0.f)&& (Frame::cosTheta(i.wo)>= 0.f)){
            v3f dif = diffuseReflectance->eval(worldData, i);
            v3f spe = specularReflectance->eval(worldData, i);
            cosA = glm::dot(i.wo, PhongBSDF::reflect(i.wi)) / ((glm::length(i.wo)) * (glm::length(PhongBSDF::reflect(i.wi))));
            if (cosA < 0){
                cosA = 0;
            }
            cos = pow(cosA, exp);
            val = (dif * INV_PI +spe * (exp+2.f) * INV_TWOPI * cos) * scale;
        }
        return val;
    }

    float pdf(const SurfaceInteraction& i) const override {
        float pdf = 0.f;
        float exp = exponent->eval(worldData,i);
        // TODO: Implement this
        v3f wr = i.frameNs.toWorld(normalize(reflect(i.wo)));
        Frame lobe(wr);
        v3f shading_world = i.frameNs.toWorld(i.wi);
        v3f dir = lobe.toLocal(shading_world);
        pdf = Warp::squareToPhongLobePdf(dir, exp);
        return pdf;
    }

    v3f sample(SurfaceInteraction& i, const v2f& _sample, float* pdf) const override {
        v3f val(0.f);
        // TODO: Implement this
        float exp = exponent->eval(worldData,i);
        v3f wr = i.frameNs.toWorld(normalize(reflect(i.wo)));
        v3f dir = Warp::squareToPhongLobe(_sample, exp); //dir is in the lobe's local frame
        Frame lobe(wr);
        v3f dir_shading_World = lobe.toWorld(dir); // sample direction toWorld
        i.wi = i.frameNs.toLocal(normalize(dir_shading_World));
        val = PhongBSDF::eval(i);
        *pdf = PhongBSDF::pdf(i);
        return val;
    }

    std::string toString() const override { return "Phong"; }
};

TR_NAMESPACE_END