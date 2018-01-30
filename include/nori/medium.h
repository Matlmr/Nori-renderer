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

#if !defined(__NORI_MEDIUM_H)
#define __NORI_MEDIUM_H

#include <nori/object.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief This is actually a isotropic homogeneous medium
 */
class Medium : public NoriObject {
public:
    Medium(const PropertyList &propList) {
        m_sigmaA = propList.getColor("absorption");
        m_sigmaS = propList.getColor("scattering");
        m_sigmaT = m_sigmaA + m_sigmaS;
        m_albedo = m_sigmaS / m_sigmaT;
        m_phaseF = 0.25 * INV_PI;
    }

    virtual std::string toString() const override {
        return tfm::format("Medium[\n"
            "]"
            );
    }

    virtual EClassType getClassType() const override { return EMedium; }
    
    Color3f transmittance(const Point3f x, const Point3f y) const;
    
    Color3f getSigmaA() const { return m_sigmaA; }
    
    Color3f getSigmaS() const { return m_sigmaS; }
    
    Color3f getSigmaT() const { return m_sigmaT; }
    
    Color3f getAlbedo() const { return m_albedo; }
    
    float getPhaseF() const { return m_phaseF; }
    
private:
    Color3f m_sigmaA, m_sigmaS, m_sigmaT, m_albedo;
    float m_phaseF;
};

NORI_NAMESPACE_END

#endif /* __NORI_MEDIUM_H */
