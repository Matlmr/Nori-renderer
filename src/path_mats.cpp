#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class PathMaterialSampling : public Integrator {
public:
    PathMaterialSampling(const PropertyList &props) {
        /* No parameters this time */
    }
    
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
      
        // Initial radiance and throughput
        Color3f Li(0.0f), t(1.f);
        Ray3f mRay = ray;
        const Emitter* env = scene->getEnvEmitter();
        Color3f BSDF(1.0f);
        
        while (true) {
            
            Intersection its;
            
            if(!scene->rayIntersect(mRay, its)) {
                // No more intersection, return current Li
                EmitterQueryRecord lRec;
                lRec.wi = mRay.d.normalized();
                // Check for envmap
                if (env != nullptr ) Li += env->eval(lRec) * BSDF;
                return Li;
            }
        
            // First check if intersection is an emitter
            if (its.mesh->isEmitter()) {
                // Create a query to get the value
                EmitterQueryRecord lRec(mRay.o, its.p, its.shFrame.n);
                Li += t*its.mesh->getEmitter()->eval(lRec);
            }
            
            // Russian roulette
            if (sampler->next1D() > std::min(t.maxCoeff(),0.999f)) {
                return Li;
            }
            t /= std::min(t.maxCoeff(),0.999f);
        
            // Use BSDF Sampling and shoot a ray in that direction
            BSDFQueryRecord bRec(its.shFrame.toLocal(-mRay.d));
            bRec.uv = its.uv;
            BSDF = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            mRay = Ray3f(its.p, its.toWorld(bRec.wo));
        
            // The sample function already returns the value divided by the pdf
            t *= BSDF;
            
        }
    }
    
    std::string toString() const {
        return "PathMaterialSampling[]";
    }
};

NORI_REGISTER_CLASS(PathMaterialSampling, "path_mats");
NORI_NAMESPACE_END
