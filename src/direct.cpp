#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class DirectIntegrator : public Integrator {
public:
    DirectIntegrator(const PropertyList &props) {
        /* No parameters this time */
    }
    
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);
        
        // Shading normal
        Normal3f n = its.shFrame.n;
        
        // Vector of lights
        std::vector<Emitter *> emitters = scene->getLights();
        
        Color3f color;
        
        // Query to get data from lights
        EmitterQueryRecord lRec;
        lRec.ref = its.p;
        Point2f sampl;
        
        for (unsigned int i = 0; i < emitters.size(); ++i) {
            
            // Call sample of emitter to fill query
            Color3f emitvalue = emitters[i]->sample(lRec, sampl);
        
            // Angle between direction from x to p and shading normal
            float cosTheta = Frame::cosTheta(its.shFrame.toLocal(lRec.wi));
        
            // Query for the BSDF
            BSDFQueryRecord bRec = BSDFQueryRecord(its.shFrame.toLocal(lRec.wi), its.shFrame.toLocal(-lRec.wi + 2*cosTheta*n), ESolidAngle);
            
            // Set the uv coordinates of the query
            bRec.uv = its.uv;
        
            // Add the color from each emitter
            color += emitvalue * (cosTheta + sqrt(pow(cosTheta,2)))/2 * its.mesh->getBSDF()->eval(bRec)*(!scene->rayIntersect(lRec.shadowRay,its));
        }
        
        // Return the sum
        return color;
    }
    
    std::string toString() const {
        return "DirectIntegrator[]";
    }
};

NORI_REGISTER_CLASS(DirectIntegrator, "direct");
NORI_NAMESPACE_END
