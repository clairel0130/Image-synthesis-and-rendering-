/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Ambient occlusion integrator
 */
struct AOIntegrator : Integrator {
    explicit AOIntegrator(const Scene& scene) : Integrator(scene) { }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);
        v3f samp(0.f);
        float cosT = 0.f;
        SurfaceInteraction i;
	    // TODO: Implement this
	    if (scene.bvh->intersect(ray, i)){
//-----------------------------------uniform spherical direction sampling-----------------------------------------------

                samp = Warp::squareToUniformSphere(sampler.next2D());
                i.wi = normalize(i.frameNs.toWorld(samp));
                cosT = Frame::cosTheta(samp);
                if (cosT >= 0.f){
                    Li = v3f(1.f * INV_PI * cosT / Warp::squareToUniformSpherePdf());
                }
                Ray shadowRay = TinyRender::Ray(i.p, i.wi);
                shadowRay.max_t = scene.aabb.getBSphere().radius / 2.f;
                if (scene.bvh->intersect(shadowRay, i)){
                    Li = v3f(0.f);
                }

//------------------------------END (Uniform Spherical)-----------------------------------------------------------------


//-----------------------------------uniform hemispherical direction sampling-------------------------------------------

//                samp = Warp::squareToUniformHemisphere(sampler.next2D());
//                i.wi = normalize(i.frameNs.toWorld(samp));
//                cosT = Frame::cosTheta(samp);
//                if (cosT >= 0.f){
//                    Li = v3f(1.f * INV_PI * cosT / Warp::squareToUniformHemispherePdf(samp));
//                }
//                Ray shadowRay = TinyRender::Ray(i.p, i.wi);
//                shadowRay.max_t = scene.aabb.getBSphere().radius / 2.f;
//                if (scene.bvh->intersect(shadowRay, i)){
//                    Li = v3f(0.f);
//                }


//------------------------------END (Uniform Hemispheical)--------------------------------------------------------------

//------------------------------Cosine-weighted hemispherical direction sampling----------------------------------------

//                samp = Warp::squareToCosineHemisphere(sampler.next2D());
//                i.wi = normalize(i.frameNs.toWorld(samp));
//                cosT = Frame::cosTheta(samp);
//                if (cosT >= 0.f){
//                    Li = v3f(1.f * INV_PI * cosT / Warp::squareToCosineHemispherePdf(samp));
//                }
//                Ray shadowRay = TinyRender::Ray(i.p, i.wi);
//                shadowRay.max_t = scene.aabb.getBSphere().radius / 2.f;
//                if (scene.bvh->intersect(shadowRay, i)){
//                    Li = v3f(0.f);
//                }


//------------------------------END (Cosine-weighted Hemispheical)------------------------------------------------------
        }
        return Li;

        //Frame frame;
    }
};

TR_NAMESPACE_END