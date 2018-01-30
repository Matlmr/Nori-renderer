/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class RoughConductor : public BSDF {
public:
    RoughConductor(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);
        /* Eta and K (default: gold) */
        std::string material = propList.getString("material","Au");
        if (material == "Au") {
            m_eta = Color3f(0.1431189557f, 0.3749570432f, 1.4424785571f);
            m_k = Color3f(3.9831604247f, 2.3857207478f, 1.6032152899f);
        } else if (material == "Cu") {
            m_eta = Color3f(0.2004376970f, 0.9240334304f, 1.1022119527f);
            m_k = Color3f(3.9129485033f, 2.4528477015f, 2.1421879552f);
        } else if (material == "Cr") {
            m_eta = Color3f(4.3696828663f, 2.9167024892f, 1.6547005413f);
            m_k = Color3f(5.2064337956f, 4.2313645277f, 3.7549467933f);
        }
    }

    /// Evaluate the microfacet normal distribution D
    float evalBeckmann(const Normal3f &m) const {
        float temp = Frame::tanTheta(m) / m_alpha,
              ct = Frame::cosTheta(m), ct2 = ct*ct;

        return std::exp(-temp*temp) 
            / (M_PI * m_alpha * m_alpha * ct2 * ct2);
    }

    /// Evaluate Smith's shadowing-masking function G1 
    float smithBeckmannG1(const Vector3f &v, const Normal3f &m) const {
        float tanTheta = Frame::tanTheta(v);

        /* Perpendicular incidence -- no shadowing/masking */
        if (tanTheta == 0.0f)
            return 1.0f;

        /* Can't see the back side from the front and vice versa */
        if (m.dot(v) * Frame::cosTheta(v) <= 0)
            return 0.0f;

        float a = 1.0f / (m_alpha * tanTheta);
        if (a >= 1.6f)
            return 1.0f;
        float a2 = a * a;

        /* Use a fast and accurate (<0.35% rel. error) rational
           approximation to the shadowing-masking function */
        return (3.535f * a + 2.181f * a2) 
             / (1.0f + 2.276f * a + 2.577f * a2);
    }

    /// Evaluate the BRDF for the given pair of directions
    virtual Color3f eval(const BSDFQueryRecord &bRec) const override {
        
        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        
        float D = evalBeckmann(wh);
        Color3f F = fresnelCond((wh.dot(bRec.wo)), m_eta, m_k);
        float G = smithBeckmannG1(bRec.wi, wh) * smithBeckmannG1(bRec.wo,wh);
        
        return D * F * G / (4.f * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo));
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    virtual float pdf(const BSDFQueryRecord &bRec) const override {
        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        if (Frame::cosTheta(bRec.wo) <= 0) {
            return 0.f;
        }
        float D = evalBeckmann(wh);
        float Jh = 1.f / (4.f * wh.dot(bRec.wo));
        return D * Frame::cosTheta(wh) * Jh;
    }

    /// Sample the BRDF
    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const override {
    
        if (Frame::cosTheta(bRec.wi) <= 0) {
            return Color3f(0.0f);
        }
        
        Point2f sample = _sample;
        Vector3f wh = Warp::squareToBeckmann(sample, m_alpha);
        bRec.wo = ((2.f * wh.dot(bRec.wi) * wh) - bRec.wi).normalized();
        
        if (Frame::cosTheta(bRec.wo) <= 0) {
            return Color3f(0.0f);
        }
        return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    virtual std::string toString() const override {
        return tfm::format(
            "RoughConductor[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "]",
            m_alpha,
            m_eta,
            m_k
        );
    }
private:
    float m_alpha;
    Color3f m_eta, m_k;
};

NORI_REGISTER_CLASS(RoughConductor, "roughconductor");
NORI_NAMESPACE_END
