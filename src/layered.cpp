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

NORI_NAMESPACE_BEGIN

class Layered : public BSDF {
public:
    Layered(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        m_thickness = propList.getFloat("thickness", 1);
        m_sigmaA = propList.getColor("absorption", Color3f(0,0,0));
    }
    
    virtual ~Layered() {
        delete m_bsdf;
    }
	
    void addChild(NoriObject *obj) override {
        switch (obj->getClassType()) {
            case EBSDF:
                m_bsdf = static_cast<BSDF *>(obj);
        }
    }
    
    virtual Color3f eval(const BSDFQueryRecord &bRec) const override {
        /*if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);*/
        if (Frame::cosTheta(bRec.wi) <= 0.f || Frame::cosTheta(bRec.wo) <= 0) {
            return Color3f(0.0f);
        }
        
        BSDFQueryRecord bRec_t(bRec);
		float through_i, through_o;
        through_i = refraction(bRec.wi, bRec_t.wi);
        through_o = refraction(bRec.wo, bRec_t.wo);
        
        Color3f BSDF = m_bsdf->eval(bRec_t);
        BSDF *= exp(-m_sigmaA * m_thickness *
            (1 / abs(Frame::cosTheta(bRec_t.wi)) + 1 / abs(Frame::cosTheta(bRec_t.wo))));

        return BSDF * (1 - through_i) * (1 - through_o);
    }

    virtual float pdf(const BSDFQueryRecord &bRec) const override {

        if (Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0) {
            return 1.0f;
        }
        
        BSDFQueryRecord bRec_t(bRec);
        float through_i = refraction(bRec.wi, bRec_t.wi);
        float through_o = refraction(bRec.wo, bRec_t.wo);

        float pdf = m_bsdf->pdf(bRec_t);
        
	    return pdf * (through_i != 1.f) * (through_o != 1.f);
    }

    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const override {
        if (Frame::cosTheta(bRec.wi) <= 0) {
            return Color3f(0.0f);
        }
        
        Vector3f transmit_i;
        float through_i = refraction(bRec.wi, transmit_i);
        
        if (sample.x() < through_i) {
            // Reflection
            bRec.wo = Vector3f(-bRec.wi.x(),-bRec.wi.y(),bRec.wi.z());
            bRec.eta = 1.0f;
            return Color3f(1.f);
        }
        else {
            // Refraction
            BSDFQueryRecord bRec_t(transmit_i);
            Color3f BSDF = m_bsdf->sample(bRec_t, sample);
            bRec = BSDFQueryRecord(bRec.wi, bRec_t.wo, EUnknownMeasure);
            bRec.eta = bRec_t.eta;
            bRec.measure = bRec_t.measure;
            Vector3f transmit_o = bRec.wo;
            
            // Compute refraction in the other direction, can't swap because const
            float cosTheta;
            float through_o = fresnelDielCos(std::abs(Frame::cosTheta(transmit_o)), cosTheta, m_intIOR, m_extIOR);
            bRec.wo = Vector3f(m_intIOR / m_extIOR * transmit_o.x(), m_intIOR / m_extIOR * transmit_o.y(), -Frame::cosTheta(transmit_o)/abs(Frame::cosTheta(transmit_o)) * cosTheta);
            
            BSDF *= exp(-m_sigmaA * m_thickness *
                    (1 / std::abs(Frame::cosTheta(transmit_i)) +
                     1 / std::abs(Frame::cosTheta(transmit_o))));
            
            return BSDF * (1 - through_i) * (1 - through_o);
        }
    }
    
    float refraction(const Vector3f &wi, Vector3f &wt) const {
        float cosTheta;
        float th = fresnelDielCos(abs(Frame::cosTheta(wi)), cosTheta, m_extIOR, m_intIOR);
        wt = Vector3f(m_extIOR / m_intIOR * wi.x(), m_extIOR / m_intIOR * wi.y(), -Frame::cosTheta(wi)/abs(Frame::cosTheta(wi)) * cosTheta);
        return th;
    }

    virtual std::string toString() const override {
        return tfm::format(
            "Layered[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "  sigma = \"%s\"\n"
            "]",
            m_intIOR, m_extIOR, m_sigmaA.toString());
    }
private:
    float m_intIOR, m_extIOR;
    float m_thickness;
    Color3f m_sigmaA;
    BSDF *m_bsdf;
};

NORI_REGISTER_CLASS(Layered, "layered");
NORI_NAMESPACE_END
