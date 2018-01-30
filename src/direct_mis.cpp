#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class DirectMultiImportanceSampling : public Integrator {
public:
    DirectMultiImportanceSampling(const PropertyList &props) {
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
        
        /* Emitter Sampling */
        
        // Random Light
        const Emitter * emitter = scene->getRandomEmitter(sampler->next1D());
        
        // Query to get data from lights
        EmitterQueryRecord lRecR_ems;
        lRecR_ems.ref = itsE.p;
        
        // Call sample of emitter to fill query
        Color3f radiance_ems = emitter->sample(lRecR_ems, sampler->next2D())*scene->getLights().size();
        float pdf_emsE = emitter->pdf(lRecR_ems);
        
        // Angle between direction from x to p and shading normal
        float cosTheta_ems = Frame::cosTheta(itsE.shFrame.toLocal(lRecR_ems.wi));
        
        // Query for the BSDF
        BSDFQueryRecord bRec_ems = BSDFQueryRecord(itsE.shFrame.toLocal(-ray.d), itsE.shFrame.toLocal(lRecR_ems.wi), ESolidAngle);
        
        // Get the BSDF value
        Color3f BSDF_ems = itsE.mesh->getBSDF()->eval(bRec_ems);
        float pdf_emsB = itsE.mesh->getBSDF()->pdf(bRec_ems);
        
        // Set the uv coordinates of the query
        bRec_ems.uv = itsE.uv;
        
        // Check the shadowray
        float obstacle = 1.f;
        Intersection its_ems;
        if (scene->rayIntersect(lRecR_ems.shadowRay,its_ems)) {
            obstacle = 0.f;
        }
        
        /* BSDF Sampling */
        
        // Use BSDF Sampling and shoot a ray in that direction
        BSDFQueryRecord bRec_mats(itsE.shFrame.toLocal(-ray.d));
        bRec_mats.uv = itsE.uv;
        Color3f BSDF_mats = itsE.mesh->getBSDF()->sample(bRec_mats, sampler->next2D());
        float pdf_matsB = itsE.mesh->getBSDF()->pdf(bRec_mats);
        Ray3f rayR = Ray3f(itsE.p, itsE.toWorld(bRec_mats.wo));
        
        // Check light emitted in that direction
        Color3f radiance_mats(0.f);
        float cosTheta_mats;
        Intersection itsR;
        float pdf_matsE = 0.f;
        if(scene->rayIntersect(rayR, itsR)) {
            if(itsR.mesh->isEmitter()) {
                EmitterQueryRecord lRecR_mats(itsE.p, itsR.p, itsR.shFrame.n);
                radiance_mats = itsR.mesh->getEmitter()->eval(lRecR_mats);
                pdf_matsE = itsR.mesh->getEmitter()->pdf(lRecR_mats);
                cosTheta_mats = Frame::cosTheta(itsE.shFrame.toLocal(lRecR_mats.wi));
            }
        }
        
        float w_em(0.f), w_mat(0.f);
        
        if ( pdf_emsE + pdf_emsB != 0.f ) {
            w_em = pdf_emsE / (pdf_emsE + pdf_emsB);
        }
        if ( pdf_matsE != 0.f + pdf_matsB != 0.f ) {
            w_mat = pdf_matsB / (pdf_matsE + pdf_matsB) ;
        }

        return Le + w_em * obstacle * radiance_ems * BSDF_ems * std::max(0.f,cosTheta_ems) + w_mat * radiance_mats * BSDF_mats;
    }
    
    std::string toString() const {
        return "DirectMultiImportanceSampling[]";
    }
};

NORI_REGISTER_CLASS(DirectMultiImportanceSampling, "direct_mis");
NORI_NAMESPACE_END
