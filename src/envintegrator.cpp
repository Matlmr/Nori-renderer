#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class EnvIntegrator : public Integrator {
public:
    EnvIntegrator(const PropertyList &props) {
        /* No parameters this time */

    }
    
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        
        EmitterQueryRecord lRec;
        lRec.wi = ray.d;

        const Emitter* emitter = scene->getRandomEmitter(sampler->next1D());
        
        return emitter->eval(lRec);
        
    }
    
    std::string toString() const {
        return "EnvIntegrator[]";
    }
};

NORI_REGISTER_CLASS(EnvIntegrator, "envintegrator");
NORI_NAMESPACE_END
