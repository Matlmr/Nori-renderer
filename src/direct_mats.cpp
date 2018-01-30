#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class DirectMaterialSampling : public Integrator {
public:
    DirectMaterialSampling(const PropertyList &props) {
        /* No parameters this time */
    }
    
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection itsE;
        if (!scene->rayIntersect(ray, itsE))
            return Color3f(0.0f);
        
        Color3f Le(0.0f);
        
        // First check if intersection is an emitter
        if (itsE.mesh->isEmitter()) {
            // Create a query to get the value
            EmitterQueryRecord lRecE(ray.o, itsE.p, itsE.shFrame.n);
            Le = itsE.mesh->getEmitter()->eval(lRecE);
        }
        
        // Use BSDF Sampling and shoot a ray in that direction
        BSDFQueryRecord bRec(itsE.shFrame.toLocal(-ray.d));
        bRec.uv = itsE.uv;
        Color3f BSDF = itsE.mesh->getBSDF()->sample(bRec, sampler->next2D());
        Ray3f rayR = Ray3f(itsE.p, itsE.toWorld(bRec.wo));
        
        // Check light emitted in that direction
        Color3f radiance(0.f);
        float cosTheta(0.f);
        Intersection itsR;
        if(scene->rayIntersect(rayR, itsR)) {
            if(itsR.mesh->isEmitter()) {
                EmitterQueryRecord lRecR(itsE.p, itsR.p, itsR.shFrame.n);
                radiance = itsR.mesh->getEmitter()->eval(lRecR);
                cosTheta = Frame::cosTheta(itsE.shFrame.toLocal(lRecR.wi));
            }
        }
        return Le + radiance*BSDF;
        
    }
    
    std::string toString() const {
        return "DirectMaterialSampling[]";
    }
};

NORI_REGISTER_CLASS(DirectMaterialSampling, "direct_mats");
NORI_NAMESPACE_END
