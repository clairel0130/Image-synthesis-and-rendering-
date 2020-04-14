/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <core/platform.h>
#include <core/integrator.h>

TR_NAMESPACE_BEGIN

/**
 * Surface normal integrator.
 */
struct NormalIntegrator : Integrator {
    explicit NormalIntegrator(const Scene& scene) : Integrator(scene) { }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f color(0.f, 1.f, 0.f);
        // Your normal integrator solution from A1 here (optional).
        SurfaceInteraction intersection;
        if (scene.bvh->intersect(ray, intersection)){
            color = glm::abs(intersection.frameNs.n);
        }   return color;
    }
};

TR_NAMESPACE_END