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

class MultiLayered : public BSDF {
public:
    MultiLayered(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }
    virtual ~MultiLayered() {
        delete m_conductor;
    }
    
    /// Add BSDF for the sub-layer (conductor)
    virtual void addChild(NoriObject *obj) override {
        switch (obj->getClassType()) {
            case EBSDF:
                m_conductor = static_cast<BSDF *>(obj);
                break;
            default:
                throw NoriException("Diffuse::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }

    /// Evaluate the microfacet normal distribution D
    float evalBeckmann(const Normal3f &m) const {
        float temp = Frame::tanTheta(m) / m_alpha,
              ct = Frame::cosTheta(m), ct2 = ct*ct;

        return std::exp(-temp*temp) 
            / (M_PI * m_alpha * m_alpha * ct2 * ct2);
    }
    
    /// Evaluate the microfacet Blinn distribution D
    float evalBlinn(const Normal3f &m) const {
        float e = m_alpha;
        return (e + 2) * 0.5f * INV_PI * pow(std::abs(Frame::cosTheta(m)),e);
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
    
    Vector3f refraction(Vector3f wi, float &R) const {
        float cosThetaI = Frame::cosTheta(wi);
        float cosThetaT;
        R = fresnelDielVec(std::abs(cosThetaI), cosThetaT, m_extIOR, m_intIOR);
        float eta1_2 = m_extIOR / m_intIOR;
        Vector3f n = Vector3f(0.0f, 0.0f, 1.0f);
        if (cosThetaI <= 0.0f) {
            eta1_2 = m_intIOR / m_extIOR;
            n *= -1;
        }
        Vector3f refracted = -eta1_2 * (wi - (wi.dot(n) * n));
        return refracted - n * sqrt(1-pow(eta1_2,2)*(1-pow(wi.dot(n),2)));
    }

    /// Evaluate the BRDF for the given pair of directions
    virtual Color3f eval(const BSDFQueryRecord &bRec) const override {
        
        /* MICROFACET MODEL */
        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        
        float D = evalBeckmann(wh);
        float F = fresnelDiel(wh.dot(bRec.wo), m_extIOR, m_intIOR);
        float G = smithBeckmannG1(bRec.wi, wh) * smithBeckmannG1(bRec.wo,wh);
        
        //return m_kd *INV_PI + m_ks * D * F * G / (4.f * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo));
        Color3f fr1 = m_kd *INV_PI + m_ks * D * F * G / (4.f * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo));
        
        /* Dielectric */
        float R1, R2;
        BSDFQueryRecord bRecR(refraction(bRec.wi, R1),refraction(bRec.wo, R2),ESolidAngle);
        Point2f sampl;
        Color3f fr2 = m_conductor->sample(bRecR,sampl) * (1 - R1) * (1 - R2);
        return fr2;
        return 0.5f * (fr1 + fr2);
        
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    virtual float pdf(const BSDFQueryRecord &bRec) const override {
        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        if (Frame::cosTheta(bRec.wo) <= 0) {
            return 0.f;
        }
        float D = evalBeckmann(wh);
        float Jh = 1.f / (4.f * wh.dot(bRec.wo));
        return m_ks * D * Frame::cosTheta(wh) * Jh + (1 - m_ks) * Frame::cosTheta(bRec.wo) * INV_PI;
    }

    /// Sample the BRDF
    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const override {
        
        if (Frame::cosTheta(bRec.wi) <= 0) {
            return Color3f(0.0f);
        }
        
        Point2f sample = _sample;
        if (sample.x() < m_ks) {
            sample.x() /= m_ks;
            Vector3f wh = Warp::squareToBeckmann(sample, m_alpha);
            bRec.wo = ((2.f * wh.dot(bRec.wi) * wh) - bRec.wi).normalized();
        } else {
            sample.x() = (sample.x() - m_ks) / (1.0f - m_ks);
            bRec.wo = Warp::squareToCosineHemisphere(sample);
        }
        
        if (Frame::cosTheta(bRec.wo) <= 0) {
            return Color3f(0.0f);
        }
        return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    virtual std::string toString() const override {
        return tfm::format(
            "MultiLayered[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
    BSDF * m_conductor;
};

NORI_REGISTER_CLASS(MultiLayered, "multilayered");
NORI_NAMESPACE_END
