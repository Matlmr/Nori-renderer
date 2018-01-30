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
        
        // Shading normal
        // Normal3f n = its.shFrame.n;
        
        // Random Light
        Point2f sample = sampler->next2D();
        const Emitter * emitter = scene->getRandomEmitter(sample.x());
        
        Color3f color;
        
        // Query to get data from lights
        EmitterQueryRecord lRec;
        lRec.ref = its.p;
        
        // Call sample of emitter to fill query
        Color3f radiance = emitter->sample(lRec, sample);
        
        // Angle between direction from x to p and shading normal
        float cosTheta = Frame::cosTheta(its.shFrame.toLocal(lRec.wi));
        
        // Query for the BSDF
        BSDFQueryRecord bRec = BSDFQueryRecord(its.shFrame.toLocal(lRec.wi), its.shFrame.toLocal(-ray.d), ESolidAngle);
        
        // Set the uv coordinates of the query
        bRec.uv = its.uv;
        
        // Add the color from each emitter
        Color3f Le = Color3f(0.f);
        if (its.mesh->isEmitter()) {
            Le = emitter->eval(lRec);
        }
        
        Color3f bsdf = Color3f(0.f);
        if (!scene->rayIntersect(lRec.shadowRay,its)) {
            bsdf = its.mesh->getBSDF()->eval(bRec);
        }
        
        color += Le + radiance * bsdf * std::max(0.f,cosTheta);
        
        // Return the sum
        return color;
    }
    
    std::string toString() const {
        return "DirectEmitterSampling[]";
    }
};

NORI_REGISTER_CLASS(DirectEmitterSampling, "direct_ems");
NORI_NAMESPACE_END
