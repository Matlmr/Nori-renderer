#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/warp.h>
#include <nori/medium.h>

NORI_NAMESPACE_BEGIN

class VolumetricPT : public Integrator {
public:
    VolumetricPT(const PropertyList &props) {
        /* No parameters this time */
    }
    
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
      
        Color3f Li(0.0f), th(1.0f);
        const Medium* medium = scene->getMedium();
        Color3f albedo = medium->getAlbedo();
        Color3f sigmaT = medium->getSigmaT();
        float phaseF = medium->getPhaseF();
        Ray3f nRay = ray;
        
        while (true) {
            
            // Find nearest surface
            Intersection its;
            float tmax;
            if (scene->rayIntersect(nRay, its)) {
                tmax = (its.p - nRay.o).norm();
            } else {
                tmax = its.t;
            }
            
            // Sample free path
            float t = -log(1 - sampler->next1D()) / sigmaT.maxCoeff();
        
            // Volume interaction
            if (t < tmax) {
                
                // Uniform sphere sampling because isotropic phase function
                Vector3f dir = Warp::squareToUniformSphere(sampler->next2D());
                nRay = Ray3f(nRay.o + t * nRay.d.normalized(), dir);
                
                // Emitter sampling
                const Emitter* emitter = scene->getRandomEmitter(sampler->next1D());
                EmitterQueryRecord lRec(nRay.o);
                Color3f Le = emitter->sample(lRec, sampler->next2D())*scene->getLights().size();
                
                // Check the shadowray
                Intersection its_sh;
                if (scene->rayIntersect(lRec.shadowRay, its_sh)) {
                    Le = 0.0f;
                }
                
                // Update throughput and add light contribution
                th *= albedo;
                Li += th * medium->transmittance(nRay.o, lRec.p) * phaseF * Le;
                
            // Surface interaction
            } else {
                
                // We hit an emitter
                if (its.mesh->isEmitter()) {
                    EmitterQueryRecord lRec(nRay.o, its.p, its.shFrame.n);
                    Color3f Le = its.mesh->getEmitter()->eval(lRec);
                    Li += th * medium->transmittance(nRay.o, its.p) * Le;

                // Normal surface
                } else {
                    
                    // Emitter sampling
                    const Emitter* emitter = scene->getRandomEmitter(sampler->next1D());
                    EmitterQueryRecord lRec(its.p);
                    Color3f Le = emitter->sample(lRec, sampler->next2D())*scene->getLights().size();
                    //float pdf_e = emitter->pdf(lRec);
                
                    BSDFQueryRecord bRec(its.shFrame.toLocal(-nRay.d), its.shFrame.toLocal(lRec.wi), ESolidAngle);
                    bRec.uv = its.uv;
                    Color3f BSDF = its.mesh->getBSDF()->eval(bRec);
                    //float pdf_b = its.mesh->getBSDF()->pdf(bRec);
                
                    // Check the shadowray
                    Intersection its_sh;
                    if (scene->rayIntersect(lRec.shadowRay, its_sh)) {
                        Le = 0.0f;
                    }
                
                    //float pdf = exp(-sigmaT.minCoeff()*tmax);
                    //float pdf = its.mesh->getBSDF()->pdf(bRec);
                    //if (pdf <= 0.f) pdf = 1.f;
                    //pdf = pdf == 0.f ? 1.f : pdf;
                    //float w = 1.f;
                    //if (pdf_e + pdf_b > 0.0f) w = pdf_e / (pdf_b + pdf_e);
                    
                    // Light contribution with evaluated BSDF
                    Li += th * BSDF * medium->transmittance(its.p, lRec.p) * Le;
                }
                
                // Update ray with direction sampled by BSDF
                BSDFQueryRecord bRec_d(its.shFrame.toLocal(-nRay.d));
                bRec_d.uv = its.uv;
                Color3f BSDF_d = its.mesh->getBSDF()->sample(bRec_d, sampler->next2D());
                nRay = Ray3f(its.p, its.toWorld(bRec_d.wo));
                
                // Throughput multiplied by BSDF of samping direction
                th *= BSDF_d;
                
            }
            // Russian roulette
            if (sampler->next1D() > std::min(th.maxCoeff(),0.99f)) {
                return Li;
            }
            th /= std::min(th.maxCoeff(),0.99f);
        }
    }
    
    std::string toString() const {
        return "VolumetricPT[]";
    }
};

NORI_REGISTER_CLASS(VolumetricPT, "vol_pt");
NORI_NAMESPACE_END
