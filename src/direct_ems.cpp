#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class DirectEmitterSampling : public Integrator {
public:
    DirectEmitterSampling(const PropertyList &props) {
        /* No parameters this time */

    }
    
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        
        /* Find the surface that is visible in the requested direction */

        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);
        
        // Random Light
        const Emitter * emitter = scene->getRandomEmitter(sampler->next1D());
        
        // Query to get data from lights
        EmitterQueryRecord lRecR;
        lRecR.ref = its.p;
        
        // Call sample of emitter to fill query
        Color3f radiance = emitter->sample(lRecR, sampler->next2D())*scene->getLights().size();
        
        // Angle between direction from x to p and shading normal
        float cosTheta = Frame::cosTheta(its.shFrame.toLocal(lRecR.wi));
        
        // Query for the BSDF
        BSDFQueryRecord bRec(its.shFrame.toLocal(-ray.d), its.shFrame.toLocal(lRecR.wi), ESolidAngle);
        // Set the uv coordinates of the query
        bRec.uv = its.uv;
        // Get the BSDF value
        Color3f BSDF = its.mesh->getBSDF()->eval(bRec);
        
        // Add the color from the emitter (first intersection)
        Color3f Le = Color3f(0.f);
        if (its.mesh->isEmitter()) {
            EmitterQueryRecord lRecE(ray.o, its.p, its.shFrame.n);
            Le = its.mesh->getEmitter()->eval(lRecE);
        }
        
        // Check the shadowray
        float obstacle = 1.f;
        if (scene->rayIntersect(lRecR.shadowRay,its)) {
            obstacle = 0.f;
        }
        
        return Le +  obstacle * radiance * BSDF * std::max(0.f,cosTheta);
    }
    
    std::string toString() const {
        return "DirectEmitterSampling[]";
    }
};

NORI_REGISTER_CLASS(DirectEmitterSampling, "direct_ems");
NORI_NAMESPACE_END
