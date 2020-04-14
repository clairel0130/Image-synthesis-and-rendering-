/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once
#include "math.h"


TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator with MIS
 */
    struct DirectIntegrator : Integrator {
        explicit DirectIntegrator(const Scene& scene) : Integrator(scene) {
            m_emitterSamples = scene.config.integratorSettings.di.emitterSamples;
            m_bsdfSamples = scene.config.integratorSettings.di.bsdfSamples;
            m_samplingStrategy = scene.config.integratorSettings.di.samplingStrategy;
        }

        static inline float balanceHeuristic(float nf, float fPdf, float ng, float gPdf) {
            float f = nf * fPdf, g = ng * gPdf;
            return f / (f + g);
        }

        void sampleSphereByCosineHemisphere(const p2f& sample,
                                            const v3f& n,
                                            const p3f& pShading,
                                            const v3f& emitterCenter,
                                            float emitterRadius,
                                            v3f& wiW,
                                            float& pdf) const {
            // TODO: Implement this

        }

        void sampleSphereByArea(const p2f& sample,
                                const p3f& pShading,
                                const v3f& emitterCenter,
                                float emitterRadius,
                                v3f& pos,
                                v3f& ne,
                                v3f& wiW,
                                float& pdf) const {
            // TODO: Implement this
            v3f dir = Warp::squareToUniformSphere(sample);
            ne = glm::normalize(dir);
            pdf = INV_FOURPI / glm::pow(emitterRadius, 2.f);
            pos = emitterCenter + emitterRadius * ne;
            wiW = pos - pShading;
        }

        void sampleSphereBySolidAngle(const p2f& sample,
                                      const p3f& pShading,
                                      const v3f& emitterCenter,
                                      float emitterRadius,
                                      v3f& wiW,
                                      float& pdf) const {
            // TODO: Implement this
            float cosThetaMax = glm::sqrt(1.f - glm::pow((emitterRadius / glm::length(emitterCenter - pShading)), 2.f));
            v3f dir = Warp::squareToUniformCone(sample, cosThetaMax);
            Frame cone(glm::normalize(emitterCenter - pShading));
            wiW = cone.toWorld(dir);
            pdf = Warp::squareToUniformConePdf(cosThetaMax);
        }

        v3f renderArea(const Ray& ray, Sampler& sampler) const {
            v3f Lr(0.f);
            v3f emCenter(0.f);
            float emRadius;
            // TODO: Implement this
            float cosT;
            float cosO;
            v3f bsdf;
            v3f pos;
            v3f ne;
            v3f wiW;
            float pdf;

            SurfaceInteraction i;
            SurfaceInteraction inter;
            // TODO: Implement this
            for(int n=0; n<m_emitterSamples; n++){
                if (scene.bvh->intersect(ray, i)){
                    if(getEmission(i) != v3f(0.f)) {
                        return getEmission(i);
                    }
                    float emPdf;
                    size_t id = selectEmitter(sampler.next(), emPdf);
                    const Emitter& em = getEmitterByID(id);
                    emCenter = scene.getShapeCenter(em.shapeID);
                    emRadius = scene.getShapeRadius(em.shapeID);
                    sampleSphereByArea(sampler.next2D(), i.p, emCenter, emRadius, pos, ne, wiW, pdf);
                    i.wi = i.frameNs.toLocal(wiW);
                    cosT = Frame::cosTheta(glm::normalize(i.wi));
                    cosO = glm::dot(glm::normalize(-wiW), glm::normalize(ne));
                    cosO = cosO / glm::pow(glm::length(wiW), 2.f);
                    if (cosO < 0){
                        cosO = 0;
                    }
                    bsdf = getBSDF(i)->eval(i);
                    Ray shadowRay = TinyRender::Ray(i.p, glm::normalize(wiW), Epsilon);

                    if (scene.bvh->intersect(shadowRay, inter)) {
                        if (getEmission(inter) != v3f(0.f)){
                            if (cosT >= 0.f) {
                                Lr += getEmission(inter)* bsdf * cosT * cosO/ pdf;
                            }
                        }
                    }
                }
            }
            Lr = Lr/m_emitterSamples;
            return Lr;
        }

        v3f renderCosineHemisphere(const Ray& ray, Sampler& sampler) const {
            v3f Lr(0.f);
            float cosT;
            v3f bsdf;
            SurfaceInteraction i;
            SurfaceInteraction inter;
            // TODO: Implement this
            for(int n=0; n<m_emitterSamples; n++){
                if (scene.bvh->intersect(ray, i)){
                    if(getEmission(i) != v3f(0.f)) {
                        return getEmission(i);
                    }

                    i.wi = Warp::squareToCosineHemisphere(sampler.next2D());
                    cosT = Frame::cosTheta(i.wi);
                    bsdf = getBSDF(i)->eval(i);
                    Ray shadowRay = TinyRender::Ray(i.p, glm::normalize(i.frameNs.toWorld(i.wi)), Epsilon);


                    if (scene.bvh->intersect(shadowRay, inter)) {
                        if (getEmission(inter) != v3f(0.f)){
                            if (cosT >= 0.f) {
                                Lr += getEmission(inter)* bsdf * cosT / Warp::squareToCosineHemispherePdf(i.wi);
                            }
                        }
                    }
                }
            }
            Lr = Lr/m_emitterSamples;
            return Lr;
        }

        v3f renderBSDF(const Ray& ray, Sampler& sampler) const {
            v3f Lr(0.f);
            float cosT;
            v3f bsdf;
            float bsdf_pdf;
            SurfaceInteraction i;
            SurfaceInteraction inter;
            // TODO: Implement this
            for(int n=0; n<m_emitterSamples; n++){
                if (scene.bvh->intersect(ray, i)){
                    if(getEmission(i) != v3f(0.f)) {
                        return getEmission(i);
                    }
                    //sample() sets BRDF,PDF,i.wi
                    bsdf = getBSDF(i)->sample(i, sampler.next2D(), &bsdf_pdf);
                    cosT = Frame::cosTheta(i.wi);
                    Ray shadowRay = TinyRender::Ray(i.p, glm::normalize(i.frameNs.toWorld(i.wi)), Epsilon);

                    if (scene.bvh->intersect(shadowRay, inter)) {
                        if (getEmission(inter) != v3f(0.f)){
                            if (cosT >= 0.f) {
                                Lr += getEmission(inter) * bsdf * cosT/ bsdf_pdf;
                            }
                        }
                    }
                }
            }
            Lr = Lr/m_emitterSamples;
            return Lr;
        }

        v3f renderSolidAngle(const Ray& ray, Sampler& sampler) const {
            v3f Lr(0.f);
            v3f emCenter(0.f);
            float emRadius;
            float cosT;
            v3f bsdf;
            v3f wiW;
            float pdf_solid;
            SurfaceInteraction i;
            SurfaceInteraction inter;
            // TODO: Implement this
            for(int n=0; n<m_emitterSamples; n++){
                if (scene.bvh->intersect(ray, i)){
                    if(getEmission(i) != v3f(0.f)) {
                        return getEmission(i);
                    }
                    float emPdf;
                    size_t id = selectEmitter(sampler.next(), emPdf);
                    const Emitter& em = getEmitterByID(id);
                    emCenter = scene.getShapeCenter(em.shapeID);
                    emRadius = scene.getShapeRadius(em.shapeID);
                    sampleSphereBySolidAngle(sampler.next2D(), i.p, emCenter, emRadius, wiW, pdf_solid);
                    i.wi = i.frameNs.toLocal(wiW);
                    cosT = Frame::cosTheta(glm::normalize(i.wi));
                    bsdf = getBSDF(i)->eval(i);
                    Ray shadowRay = TinyRender::Ray(i.p, glm::normalize(wiW), Epsilon);

                    if (scene.bvh->intersect(shadowRay, inter)) {
                        if (getEmission(inter) != v3f(0.f)){
                            if (cosT >= 0.f) {
                                Lr += getEmission(inter)* bsdf * cosT / pdf_solid/ emPdf;
                            }
                        }
                    }
                }
            }
            Lr = Lr/m_emitterSamples;
            return Lr;
        }

        v3f renderMIS(const Ray& ray, Sampler& sampler) const {
            v3f Lr_S(0.f);
            v3f Lr_B(0.f);
            v3f Lr (0.f);
            v3f emCenter(0.f);
            float emRadius;
            float cosT;
            v3f bsdf;
            v3f wiW;
            float pdf_solid;
            float bsdf_pdf;
            float balance;
            SurfaceInteraction i;
            SurfaceInteraction inter;
            // TODO: Implement this
//-----------------------------------------Solid angle-----------------------------------------------------
            for(int n=0; n<m_emitterSamples; n++){
                if (scene.bvh->intersect(ray, i)){
                    if(getEmission(i) != v3f(0.f)) {
                        return getEmission(i);
                    }
                    float emPdf;
                    size_t id = selectEmitter(sampler.next(), emPdf);
                    const Emitter& em = getEmitterByID(id);
                    emCenter = scene.getShapeCenter(em.shapeID);
                    emRadius = scene.getShapeRadius(em.shapeID);
                    sampleSphereBySolidAngle(sampler.next2D(), i.p, emCenter, emRadius, wiW, pdf_solid);
                    i.wi = i.frameNs.toLocal(wiW);
                    cosT = Frame::cosTheta(glm::normalize(i.wi));
                    bsdf = getBSDF(i)->eval(i);
                    Ray shadowRay = TinyRender::Ray(i.p, glm::normalize(wiW), Epsilon);
                    balance = balanceHeuristic(m_emitterSamples, emPdf * pdf_solid , m_bsdfSamples, getBSDF(i)->pdf(i));
                    if (scene.bvh->intersect(shadowRay, inter)) {
                        if (getEmission(inter) != v3f(0.f)){
                            if (cosT >= 0.f) {
                                Lr_S += getEmission(inter)* bsdf * cosT * balance / pdf_solid/ emPdf;
                            }
                        }
                    }
                }
            }
            if (m_emitterSamples != 0){
                Lr_S = Lr_S /m_emitterSamples;
            }
//-----------------------------------------END-----------------------------------------------

//-----------------------------------------BSDF-----------------------------------------------------


            for(int n=0; n<m_bsdfSamples; n++){
                if (scene.bvh->intersect(ray, i)){
                    if(getEmission(i) != v3f(0.f)) {
                        return getEmission(i);
                    }

                    //sample() sets BRDF,PDF,i.wi
                    bsdf = getBSDF(i)->sample(i, sampler.next2D(), &bsdf_pdf);
                    cosT = Frame::cosTheta(i.wi);
                    Ray shadowRay = TinyRender::Ray(i.p, glm::normalize(i.frameNs.toWorld(i.wi)), Epsilon);

                    if (scene.bvh->intersect(shadowRay, inter)) {
                        if (getEmission(inter) != v3f(0.f)){
                            if (cosT >= 0.f) {
                                float emPdf;
                                emPdf = 1.f/scene.emitters.size();
                                const Emitter& em = getEmitterByID(getEmitterIDByShapeID(inter.shapeID));
                                emCenter = scene.getShapeCenter(em.shapeID);
                                emRadius = scene.getShapeRadius(em.shapeID);
                                sampleSphereBySolidAngle(sampler.next2D(), i.p, emCenter, emRadius, wiW, pdf_solid);
                                float balance_B = balanceHeuristic(m_bsdfSamples, bsdf_pdf, m_emitterSamples, pdf_solid * emPdf);
                                Lr_B += getEmission(inter) * bsdf * cosT * balance_B / bsdf_pdf;
                            }
                        }
                    }
                }
            }
            if (m_bsdfSamples!= 0){
                Lr_B = Lr_B/m_bsdfSamples;
            }


//-----------------------------------------END-----------------------------------------------
            Lr = Lr_B + Lr_S;
            return Lr;
        }

        v3f render(const Ray& ray, Sampler& sampler) const override {
            if (m_samplingStrategy == "mis")
                return this->renderMIS(ray, sampler);
            else if (m_samplingStrategy == "area")
                return this->renderArea(ray, sampler);
            else if (m_samplingStrategy == "solidAngle")
                return this->renderSolidAngle(ray, sampler);
            else if (m_samplingStrategy == "cosineHemisphere")
                return this->renderCosineHemisphere(ray, sampler);
            else if (m_samplingStrategy == "bsdf")
                return this->renderBSDF(ray, sampler);
            std::cout << "Error: wrong strategy" << std::endl;
            exit(EXIT_FAILURE);
        }

        size_t m_emitterSamples;     // Number of emitter samples
        size_t m_bsdfSamples;        // Number of BSDF samples
        string m_samplingStrategy;   // Sampling strategy to use
    };

TR_NAMESPACE_END
