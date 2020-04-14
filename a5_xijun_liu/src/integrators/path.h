/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advankced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Path tracer integrator
 */
    struct PathTracerIntegrator : Integrator {
        explicit PathTracerIntegrator(const Scene& scene) : Integrator(scene) {
            m_isExplicit = scene.config.integratorSettings.pt.isExplicit;
            m_maxDepth = scene.config.integratorSettings.pt.maxDepth;
            m_rrDepth = scene.config.integratorSettings.pt.rrDepth;
            m_rrProb = scene.config.integratorSettings.pt.rrProb;
        }


        v3f renderImplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit, int time) const {
            v3f Li(0.f);
            float cosT;
            SurfaceInteraction inter;
            v3f bsdf(0.f);
            float bsdf_pdf;

            // TODO: Implement this
            if (scene.bvh -> intersect(ray, hit)) {
                float back = glm::dot(v3f(0.f, 0.f, 1.f), hit.wo);
                if ((getEmission(hit) != v3f(0.f)) && back > 0.f) {
                    return getEmission(hit);
                }
                bsdf = getBSDF(hit)->sample(hit, sampler.next2D(), &bsdf_pdf);
                cosT = Frame::cosTheta(hit.wi);
                Ray shadowRay = TinyRender::Ray(hit.p, glm::normalize(hit.frameNs.toWorld(hit.wi)), Epsilon);

                if (scene.bvh->intersect(shadowRay, inter) && (time < m_maxDepth)) {
                    time = time + 1;
                    if (cosT >= 0.f) {
                        if (bsdf_pdf > 0) {
                            Li = bsdf * cosT * renderImplicit(shadowRay, sampler, inter, time) / bsdf_pdf;
                        }
                    }
                }
            }
            return Li;
        }

        v3f renderExplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit, int time) const {
            v3f Li(0.f);
            float RR = 1.f;
            // TODO: Implement this
            float back = glm::dot(v3f(0.f, 0.f, 1.f), hit.wo);
            if ((getEmission(hit) != v3f(0.f)) && back > 0.f){
                return getEmission(hit);
            }

            if(time < m_maxDepth || m_maxDepth == -1){
                time++;

                if (time >= m_rrDepth) {
                    if(sampler.next() > m_rrProb) {
                        return v3f(0.f);
                    }
                    else {
                        RR = m_rrProb;
                    }
                }
                //direct light
                Li += direct(hit,sampler);

                //do indirect light
                if (m_maxDepth>1 ||  m_maxDepth == -1) {
                    Li += indirect(hit,sampler,time)/RR;
                }
            }
            return Li;
        }


        v3f direct(SurfaceInteraction& hit, Sampler& sampler) const {
            v3f direct(0.f);

            float emPdf;
            size_t id = selectEmitter(sampler.next(), emPdf);
            const Emitter& em = getEmitterByID(id);
            v3f n_em;
            v3f pos_em;
            float pdf_em;

            float cosT;
            float cosO;
            v3f bsdf;
            SurfaceInteraction inter;

            sampleEmitterPosition(sampler, em, n_em, pos_em, pdf_em);
            v3f wiW = pos_em - hit.p;
            hit.wi = hit.frameNs.toLocal(glm::normalize(wiW));
            cosT = Frame::cosTheta(glm::normalize(hit.wi));
            cosO = glm::dot(glm::normalize(-wiW), glm::normalize(n_em));
            cosO = cosO / glm::pow(glm::length(wiW), 2.f);
            if (cosO < 0){
                cosO = 0;
            }
            bsdf = getBSDF(hit)->eval(hit);
            Ray shadowRay = TinyRender::Ray(hit.p, glm::normalize(wiW), Epsilon);

            if (scene.bvh->intersect(shadowRay, inter)) {
                if (getEmission(inter) != v3f(0.f)){
                    if (cosT >= 0.f) {
                        direct = getEmission(inter)* bsdf * cosT * cosO/ pdf_em;
                    }
                }
            }
            return direct;
        }


        v3f indirect(SurfaceInteraction& hit, Sampler& sampler, int time) const {
            v3f indirect_L(0.f);
            float cosT;
            SurfaceInteraction inter;
            v3f bsdf(0.f);
            float bsdf_pdf;

            bsdf = getBSDF(hit)->sample(hit, sampler.next2D(), &bsdf_pdf);
            cosT = Frame::cosTheta(glm::normalize(hit.wi));
            Ray shadowRay = TinyRender::Ray(hit.p, glm::normalize(hit.frameNs.toWorld(hit.wi)), Epsilon);

            if (scene.bvh->intersect(shadowRay, inter)) {
                if(getEmission(inter) != v3f(0.f)){
                    return v3f(0.f);
                }
                if (cosT > 0.f) {
                    if (bsdf_pdf > 0) {
                        indirect_L = bsdf * cosT * renderExplicit(shadowRay, sampler, inter, time) / bsdf_pdf;
                    }
                }
            }

            return indirect_L;
        }

        v3f render(const Ray& ray, Sampler& sampler) const override {
            Ray r = ray;
            SurfaceInteraction hit;

            if (scene.bvh->intersect(r, hit)) {
                if (m_isExplicit)
                    return this->renderExplicit(ray, sampler, hit, 0);
                else
                    return this->renderImplicit(ray, sampler, hit, 0);
            }
            return v3f(0.0);
        }

        int m_maxDepth;     // Maximum number of bounces
        int m_rrDepth;      // When to start Russian roulette
        float m_rrProb;     // Russian roulette probability
        bool m_isExplicit;  // Implicit or explicit
    };

TR_NAMESPACE_END
