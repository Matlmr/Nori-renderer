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

/// Ideal conductor BSDF
class Conductor : public BSDF {
public:
    Conductor(const PropertyList &propList) {
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

    virtual Color3f eval(const BSDFQueryRecord &) const override {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    virtual float pdf(const BSDFQueryRecord &) const override {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &) const override {
        //std::cout << "how is good\n";
        float cosTheta = Frame::cosTheta(bRec.wi);
        if (cosTheta < 0.0f) {
            return Color3f(0.0f);
        }
        bRec.measure = EDiscrete;
        bRec.wo = Vector3f(-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z());
        bRec.eta = 1.0f;
        return fresnelCond(cosTheta,m_eta,m_k);
    }

    virtual std::string toString() const override {
        return tfm::format(
            "Conductor[\n"
            "  m_eta = %f,\n"
            "  m_k = %f\n"
            "]",
            m_eta, m_k);
    }
private:
    Color3f m_eta, m_k;
};

NORI_REGISTER_CLASS(Conductor, "conductor");
NORI_NAMESPACE_END
