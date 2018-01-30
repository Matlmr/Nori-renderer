#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class PathMultiImportanceSampling : public Integrator {
public:
    PathMultiImportanceSampling(const PropertyList &props) {
        /* No parameters this time */
    }
    
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        
        // Initial radiance and throughput
        Color3f Li(0.0f), t(1.f);
        Ray3f mRay = ray;
        
        const Emitter* env = scene->getEnvEmitter();
        Color3f BSDF(1.0f);
        
        float pdf_matsB(0.f), pdf_matsE(0.f);
        float pdf_emsB(0.f), pdf_emsE(0.f);
        float cosTheta_ems;
        float w_mat(1.0f), w_em(1.0f);
        
        while (true) {
            Intersection its;
            if(!scene->rayIntersect(mRay, its)) {
                // No more intersection, return current Li
                EmitterQueryRecord lRec;
                lRec.wi = mRay.d.normalized();
                if (env != nullptr) Li += w_mat * env->eval(lRec) * t;
                return Li;
            }
            /*             *\
            Material Sampling
            \*             */
            // First check if intersection is an emitter
            Color3f radiance_mats = 0.f;
            if (its.mesh->isEmitter()) {
                // Create a query to get the value
                EmitterQueryRecord lRec_mats(mRay.o, its.p, its.shFrame.n);
                radiance_mats = t*its.mesh->getEmitter()->eval(lRec_mats);
            }
            Li += w_mat * radiance_mats;
            
            // Russian roulette
            if (sampler->next1D() > std::min(t.maxCoeff(),0.99f)) {
                return Li;
            }
            t /= std::min(t.maxCoeff(),0.99f);
            
            /*            *\
            Emitter Sampling
            \*            */
            // Random Light
            const Emitter* emitter = scene->getRandomEmitter(sampler->next1D());
            // Query to get data from lights
            EmitterQueryRecord lRec_ems(its.p);
            // Call sample of emitter to fill query
            Color3f radiance_ems = emitter->sample(lRec_ems, sampler->next2D())*scene->getLights().size();
            pdf_emsE = emitter->pdf(lRec_ems);
            // Angle between direction from x to p and normal
            cosTheta_ems = Frame::cosTheta(its.shFrame.toLocal(lRec_ems.wi));
            // Query for the BSDF
            BSDFQueryRecord bRec_ems(its.shFrame.toLocal(-mRay.d),its.shFrame.toLocal(lRec_ems.wi), ESolidAngle);
            bRec_ems.uv = its.uv;
            // Get the BSDF value
            Color3f BSDF_ems = its.mesh->getBSDF()->eval(bRec_ems);
            pdf_emsB = its.mesh->getBSDF()->pdf(bRec_ems);
            
            // Check the shadow ray
            Intersection its_ems;
            if (scene->rayIntersect(lRec_ems.shadowRay, its_ems)) {
                radiance_ems = 0.0f;
            }
            if (pdf_emsE + pdf_emsB != 0.0f) {
                w_em = pdf_emsE / (pdf_emsE + pdf_emsB);
            }

            Li += w_em * t * radiance_ems * BSDF_ems * std::max(0.f, cosTheta_ems);
            
            // Use BSDF Sampling and shoot a ray in that direction
            BSDFQueryRecord bRec(its.shFrame.toLocal(-mRay.d));
            bRec.uv = its.uv;
            Color3f BSDF = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            mRay = Ray3f(its.p, its.toWorld(bRec.wo));
            pdf_matsB = its.mesh->getBSDF()->pdf(bRec);
            // The sample function already returns the value divided by the pdf
            t *= BSDF;
            
            Intersection itsR;
            if (scene->rayIntersect(mRay, itsR)) {
                if (itsR.mesh->isEmitter()) {
                    EmitterQueryRecord lRec_R(its.p, itsR.p, itsR.shFrame.n);
                    pdf_matsE = itsR.mesh->getEmitter()->pdf(lRec_R);
                    if (pdf_matsE + pdf_matsB != 0.f) {
                        w_mat = pdf_matsB / (pdf_matsB + pdf_matsE);
                    }
                }
            } else if (env != nullptr) {
                EmitterQueryRecord lRec_R;
                lRec_R.wi = mRay.d.normalized();
                pdf_matsE = env->pdf(lRec_R);
                if (pdf_matsE + pdf_matsB != 0.f) {
                    w_mat = pdf_matsB / (pdf_matsB + pdf_matsE);
                }
            }
                
            if (bRec.measure == EDiscrete) {
                w_mat = 1.f;
                w_em = 0.f;
            }
        }
    }
    std::string toString() const {
        return "PathMultiImportanceSampling[]";
    }
};

NORI_REGISTER_CLASS(PathMultiImportanceSampling, "path_mis");
NORI_NAMESPACE_END
