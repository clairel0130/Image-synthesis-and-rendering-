/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Simple direct illumination integrator.
 */
struct SimpleIntegrator : Integrator {
    explicit SimpleIntegrator(const Scene& scene) : Integrator(scene) { }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);
	    // TODO: Add previous assignment code (if needed)
        v3f position = scene.getFirstLightPosition();
        v3f intensity = scene.getFirstLightIntensity();
        v3f disV(0.f);
        float dist = 0.f;
        SurfaceInteraction i;

        if ( scene.bvh->intersect(ray, i)){
            Ray shadowRay = TinyRender::Ray(i.p, normalize((position - i.p)));
            shadowRay.max_t = distance(i.p, position);
            if (!scene.bvh->intersect(shadowRay, i)){
                disV = position - i.p;
                dist = disV.x*disV.x + disV.y*disV.y + disV.z*disV.z;
                i.wi = normalize(i.frameNs.toLocal (position - i.p));
                Li = (intensity/dist) * (getBSDF(i)->eval(i));
            }

        }
        return Li;
    }
};

TR_NAMESPACE_END