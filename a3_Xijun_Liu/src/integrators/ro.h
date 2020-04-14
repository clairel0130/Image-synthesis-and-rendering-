/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once
#include <random>

TR_NAMESPACE_BEGIN

/**
 * Reflection occlusion integrator
 */
struct ROIntegrator : Integrator {

    float m_exponent;

    explicit ROIntegrator(const Scene& scene) : Integrator(scene) {
        m_exponent = scene.config.integratorSettings.ro.exponent;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);
        v3f wr(0.f);
        v3f wr_world(0.f);
        v3f sample(0.f);
        v3f sample_world(0.f);
        v3f sample_shadingPoint(0.f);
        float cosAlpha = 0.f;
        float cos_i = 0.f;
        float f_r = 0.f;
	    // TODO: Implement this
        SurfaceInteraction i;
        if(scene.bvh->intersect(ray, i)){
            wr = reflect(i.wo);
            wr_world = i.frameNs.toWorld(wr);
            // create a Phong lobe frame
            Frame lobe(wr_world);
            // sample a point on the lobe, convert this sample direction in the world and local frame coordinate
            sample = Warp::squareToPhongLobe(sampler.next2D(), m_exponent);
            sample_world = normalize(lobe.toWorld(sample));
            sample_shadingPoint = normalize(i.frameNs.toLocal(sample_world));
            cos_i = Frame::cosTheta(sample_shadingPoint);

            cosAlpha = dot(normalize(wr), sample_shadingPoint);
            if (cosAlpha < 0){
               cosAlpha = 0;
            }
            f_r = (m_exponent + 2.f) * INV_TWOPI * glm::pow(cosAlpha, m_exponent);
            if (Frame::cosTheta(i.wo) >= 0.f && cos_i >= 0.f){
                Li = v3f(f_r * cos_i / Warp::squareToPhongLobePdf(sample, m_exponent));
            }
            Ray shadowRay = TinyRender::Ray(i.p, sample_world, Epsilon);
            if (scene.bvh->intersect(shadowRay, i)) {
                Li = v3f(0.f);
            }
        }
        return Li;
    }
};

TR_NAMESPACE_END