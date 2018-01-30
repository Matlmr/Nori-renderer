#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class PointLight : public Emitter {
public:
    PointLight(const PropertyList & propList) {
        m_lightPos = propList.getPoint3("position", Point3f());
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
        
        lRec.shadowRay = Ray3f(lRec.ref, lRec.wi, Epsilon,(lRec.p-lRec.ref).norm());
        
        return m_power / (4 * M_PI * (m_lightPos-lRec.ref).squaredNorm() );
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
        return m_power;
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
        return 1.f;
    }
    
    virtual std::string toString() const {
        return "PointLight[]";
    }
    
    Point3f getPosition() const {
        return m_lightPos;
    }

    
private:
    Point3f m_lightPos;
    Color3f m_power;

};

NORI_REGISTER_CLASS(PointLight, "point");
NORI_NAMESPACE_END
