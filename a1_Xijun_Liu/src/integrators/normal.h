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
        v3f color = v3f (0.f, 0.f, 0.f);
	// TODO: Implement this (5 ish)
	// 1. scene bvh boolean intersect (whether ray and surface interaction)
	// 2. true : normal vector into color(integrator->rgb.data)
        SurfaceInteraction intersection;
	    if (scene.bvh->intersect(ray, intersection)){
	        color = glm::abs(intersection.frameNs.n);
    	}   return color;
    }
};

TR_NAMESPACE_END