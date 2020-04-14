/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <GL/glew.h>
#include <functional>
#include "platform.h"
#include "core.h"
#include "utils.h"
#include "cpptoml.h"
#include "tiny_obj_loader.h"
#include "camera.h"

TR_NAMESPACE_BEGIN

    inline float safeSqrt(float v){
        return std::sqrt(std::max(float(0), v));
    }

/**
 * Computes barycentric coordinates.
 */
    template<class T>
    inline T barycentric(const T& a, const T& b, const T& c, const float u, const float v) {
        return a * (1 - u - v) + b * u + c * v;
    }

/**
 * Restricts a value to a given interval.
 */
    template<class T>
    inline T clamp(T v, T min, T max) {
        return std::min(std::max(v, min), max);
    }

/**
 * Checks if vector is zero.
 */
    inline bool isZero(const v3f v) {
        return glm::dot(v, v) < Epsilon;
    }

/**
 * Generates coordinate system.
 */
    inline void coordinateSystem(const v3f& a, v3f& b, v3f& c) {
        if (std::abs(a.x) > std::abs(a.y)) {
            float invLen = 1.f / std::sqrt(a.x * a.x + a.z * a.z);
            c = v3f(a.z * invLen, 0.f, -a.x * invLen);
        } else {
            float invLen = 1.f / std::sqrt(a.y * a.y + a.z * a.z);
            c = v3f(0.f, a.z * invLen, -a.y * invLen);
        }
        b = glm::cross(c, a);
    }

/**
 * Converts RGB value to luminance.
 */
    inline float getLuminance(const v3f& rgb) {
        return glm::dot(rgb, v3f(0.212671f, 0.715160f, 0.072169f));
    }

/**
 * Pseudo-random sampler (Mersenne Twister 19937) structure.
 */
    struct Sampler {
        std::mt19937 g;
        std::uniform_real_distribution<float> d;
        explicit Sampler(int seed) {
            g = std::mt19937(seed);
            d = std::uniform_real_distribution<float>(0.f, 1.f);
        }
        float next() { return d(g); }
        p2f next2D() { return {d(g), d(g)}; }
        void setSeed(int seed) {
            g.seed(seed);
            d.reset();
        }
    };

/**
 * 1D discrete distribution.
 */
    struct Distribution1D {
        std::vector<float> cdf{0};
        bool isNormalized = false;

        inline void add(float pdfVal) {
            cdf.push_back(cdf.back() + pdfVal);
        }

        size_t size() {
            return cdf.size() - 1;
        }

        float normalize() {
            float sum = cdf.back();
            for (float& v : cdf) {
                v /= sum;
            }
            isNormalized = true;
            return sum;
        }

        inline float pdf(size_t i) const {
            assert(isNormalized);
            return cdf[i + 1] - cdf[i];
        }

        int sample(float sample) const {
            assert(isNormalized);
            const auto it = std::upper_bound(cdf.begin(), cdf.end(), sample);
            return clamp(int(distance(cdf.begin(), it)) - 1, 0, int(cdf.size()) - 2);
        }
    };


/**
 * Warping functions.
 */
    namespace Warp {

        inline v2f squareToUniformTriangle(const p2f& sample) {
            v2f v(0.f);
            float u = std::sqrt(1.f - sample.x);
            v = {1 - u, u * sample.y};
            return v;
        }

        inline v3f squareToUniformSphere(const p2f& sample) {
            v3f v(0.f);
            float r = 0;
            float phi = 0;
            // TODO: Implement this
            v.z = 2 * sample.x - 1;
            r = sqrt(1 - pow(v.z, 2.f));
            phi = 2 * M_PI * sample.y;
            v.x = r * cos(phi);
            v.y = r * sin(phi);
            return v;

        }

        inline float squareToUniformSpherePdf() {
            float pdf = 0.f;
            // TODO: Implement this
            pdf = INV_FOURPI;
            return pdf;
        }

        inline v3f squareToUniformHemisphere(const p2f& sample) {
            v3f v(0.f);
            float r = 0;
            float phi = 0;
            // TODO: Implement this
            v.z = glm::abs(2 * sample.x - 1);
            r = sqrt(1 - pow(v.z, 2.f));
            phi = 2 * M_PI * sample.y;
            v.x = r * cos(phi);
            v.y = r * sin(phi);
            return v;
        }

        inline float squareToUniformHemispherePdf(const v3f& v) {
            float pdf = 0.f;
            // TODO: Implement this
            pdf = INV_TWOPI;
            return pdf;
        }

        inline v2f squareToUniformDiskConcentric(const p2f& sample) {
            v2f v(0.f);
            // TODO: Implement this (optional)
            float r = sample.x;
            float theta = 2 * M_PI * sample.y;
            // convert from polar coordinate to Cartesian coordinate
            v.x = glm::sqrt(r) * cos(theta);
            v.y = glm::sqrt(r) * sin(theta);
            return v;
        }

        inline v3f squareToCosineHemisphere(const p2f& sample) {
            v3f v(0.f);
            // TODO: Implement this
            v2f disk = Warp::squareToUniformDiskConcentric(sample);
            v.x = disk.x;
            v.y = disk.y;
            v.z = sqrt(1 - pow(v.x, 2.f) - pow(v.y, 2.f));
            return v;
        }

        inline float squareToCosineHemispherePdf(const v3f& v) {
            float pdf = 0.f;
            // TODO: Implement this
            pdf = v.z * INV_PI;
            return pdf;
        }

        inline v3f squareToPhongLobe(const p2f& sample, float exponent) {
            v3f v(0.f);
            // TODO: Implement this
            float theta = glm::acos(pow((1.f - sample.x), 1.f/(exponent + 2.f)));
            float phi = 2.f * M_PI * sample.y;
            v.x = glm::sin(theta) * glm::cos(phi);
            v.y = glm::sin(theta) * glm::sin(phi);
            v.z = glm::cos(theta);
            return v;
        }

        inline float squareToPhongLobePdf(const v3f& v, float exponent) {
            float pdf = 0.f;
            // TODO: Implement this
            if(v.z <= 0.f){
                return 0.f;
            }
            pdf = (exponent + 2.f) * INV_TWOPI * pow(v.z, exponent);
            return pdf;
        }

        inline v3f squareToUniformCone(const p2f& sample, float cosThetaMax) {
            v3f v(0.f);
            // TODO: Implement this
            float cosTheta = (1.f - sample.x) + sample.x * cosThetaMax;
            float sinTheta = glm::sqrt(1.f - glm::pow(cosTheta, 2.f));
            float phi = sample.y * M_PI * 2.f;
            v.x = glm::cos(phi) * sinTheta;
            v.y = glm::sin(phi) * sinTheta;
            v.z = cosTheta;
            return v;
        }

        inline float squareToUniformConePdf(float cosThetaMax) {
            float pdf = 0.f;
            // TODO: Implement this
            pdf = INV_TWOPI / (1.f - cosThetaMax);
            return pdf;
        }

    }

TR_NAMESPACE_END
