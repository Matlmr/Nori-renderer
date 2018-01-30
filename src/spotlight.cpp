#include <nori/emitter.h>
#include <nori/warp.h>
#include <nori/scene.h>
#include <nori/integrator.h>

NORI_NAMESPACE_BEGIN

class SpotLight : public Emitter {
public:
    SpotLight(const PropertyList & propList) {
        m_lightPos = propList.getPoint3("position", Point3f());
        Point3f target = propList.getPoint3("target", Point3f());
        m_lightDir = (target-m_lightPos).normalized();
        m_cosCone = cos(M_PI/180*propList.getFloat("cone"));
        m_cosFalloff = cos(M_PI/180*propList.getFloat("falloff"));
        m_power = propList.getColor("power", Color3f());
    }
    
    /**
     * \brief Sample the emitter and return the importance weight (i.e. the
     * value of the Emitter divided by the probability density
     * of the sample with respect to solid angles).
     *
     * \param lRec    An emitter query record (only ref is needed)
     * \param sample  A uniformly distributed sample on \f$[0,1]^2\f$
     *
     * \return The emitter value divided by the probability density of the sample.
     *         A zero value means that sampling failed.
     */
    virtual Color3f sample(EmitterQueryRecord &lRec, const Point2f &sample) const {
        lRec.p = m_lightPos;
        lRec.wi = (lRec.p - lRec.ref).normalized();
        lRec.n = m_lightDir.normalized();
        lRec.shadowRay = Ray3f(lRec.ref, lRec.wi, Epsilon, (lRec.p - lRec.ref).norm());
        
        return eval(lRec) / pdf(lRec);
    }
    
    float falloff(Vector3f vec) const {
        float cosTheta = vec.normalized().dot(m_lightDir);
        if (cosTheta < m_cosCone) return 0.f;
        if (cosTheta > m_cosFalloff) return 1.f;
        float delta = (cosTheta - m_cosCone) / (m_cosFalloff - m_cosCone);
        return pow(delta,4);
    }
    
    /**
     * \brief Evaluate the emitter
     *
     * \param lRec
     *     A record with detailed information on the emitter query
     * \return
     *     The emitter value, evaluated for each color channel
     */
    virtual Color3f eval(const EmitterQueryRecord &lRec) const {
        return m_power * falloff(-lRec.wi) / (2 * M_PI * (1 - 0.5 * (m_cosFalloff + m_cosCone)));
    }
    
    /**
     * \brief Compute the probability of sampling \c lRec.p.
     *
     * This method provides access to the probability density that
     * is realized by the \ref sample() method.
     *
     * \param lRec
     *     A record with detailed information on the emitter query
     *
     * \return
     *     A probability/density value
     */
    virtual float pdf(const EmitterQueryRecord &lRec) const{
        float cosTheta = abs((-lRec.wi).dot(m_lightDir));
        return Warp::squareToUniformSphereCapPdf(Vector3f(0,0,1),m_cosCone) * (lRec.p - lRec.ref).squaredNorm() / cosTheta;
    }
        
    virtual std::string toString() const {
        return "SpotLight[]";
    }
    
    Point3f getPosition() const {
        return m_lightPos;
    }

    
private:
    Point3f m_lightPos;
    Color3f m_power;
    Point3f m_lightDir;
    float m_cosCone;
    float m_cosFalloff;

};

NORI_REGISTER_CLASS(SpotLight, "spot");
NORI_NAMESPACE_END

/*
 #include <nori/emitter.h>
 #include <nori/warp.h>
 
 NORI_NAMESPACE_BEGIN
 
 class SpotLight : public Emitter {
 public:
 SpotLight(const PropertyList & propList) {
 m_lightPos = propList.getPoint3("position", Point3f());
 Point3f dest = propList.getPoint3("target", Point3f());
 m_lightDir = (dest-m_lightPos).normalized();
 m_cosCone = cos(M_PI/180*propList.getFloat("cone"));
 m_cosFalloff = cos(M_PI/180*propList.getFloat("falloff"));
 m_power = propList.getColor("power", Color3f());
 }
 

virtual Color3f sample(EmitterQueryRecord &lRec, const Point2f &sample) const {
    lRec.p = m_lightPos;
    lRec.wi = (lRec.p - lRec.ref).normalized();
    lRec.n = m_lightDir.normalized();
    lRec.shadowRay = Ray3f(lRec.ref, lRec.wi, Epsilon, lRec.wi.norm());
    
    return eval(lRec) / pdf(lRec);
}

float falloff(Vector3f vec) const {
    float cosTheta = vec.dot(m_lightDir);
    if (cosTheta < m_cosCone) return 0.f;
    if (cosTheta > m_cosFalloff) return 1.f;
    float delta = (cosTheta - m_cosCone) / (m_cosFalloff - m_cosCone);
    return pow(delta,4);
}


virtual Color3f eval(const EmitterQueryRecord &lRec) const {
    return m_power * falloff(-lRec.wi) / (2 * M_PI * (1 - 0.5 * (m_cosFalloff + m_cosCone)));
}


virtual float pdf(const EmitterQueryRecord &lRec) const{
    float cosTheta = (-lRec.wi).dot(m_lightDir);
    return Warp::squareToUniformSphereCapPdf(Vector3f(0,0,1),m_cosCone) * (lRec.p - lRec.ref).squaredNorm() / cosTheta;
}

virtual std::string toString() const {
    return "SpotLight[]";
}

Point3f getPosition() const {
    return m_lightPos;
}


private:
Point3f m_lightPos;
Color3f m_power;
Point3f m_lightDir;
float m_cosCone;
float m_cosFalloff;

};

NORI_REGISTER_CLASS(SpotLight, "spot");
NORI_NAMESPACE_END
*/
