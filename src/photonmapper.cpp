/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/scene.h>
#include <nori/photon.h>

NORI_NAMESPACE_BEGIN

class PhotonMapper : public Integrator {
public:
    /// Photon map data structure
    typedef PointKDTree<Photon> PhotonMap;

    PhotonMapper(const PropertyList &props) {
        /* Lookup parameters */
        m_photonCount  = props.getInteger("photonCount", 1000000);
        m_photonRadius = props.getFloat("photonRadius", 0.0f /* Default: automatic */);
    }

    virtual void preprocess(const Scene *scene) override {
        cout << "Gathering " << m_photonCount << " photons .. ";
        cout.flush();

        /* Create a sample generator for the preprocess step */
        Sampler *sampler = static_cast<Sampler *>(
            NoriObjectFactory::createInstance("independent", PropertyList()));

        /* Allocate memory for the photon map */
        m_photonMap = std::unique_ptr<PhotonMap>(new PhotonMap());
        m_photonMap->reserve(m_photonCount);

		/* Estimate a default photon radius */
		if (m_photonRadius == 0)
			m_photonRadius = scene->getBoundingBox().getExtents().norm() / 500.0f;

	

		/* How to add a photon?
		 * m_photonMap->push_back(Photon(
		 *	Point3f(0, 0, 0),  // Position
		 *	Vector3f(0, 0, 1), // Direction
		 *	Color3f(1, 2, 3)   // Power
		 * ));
		 */

		// put your code to trace photons here
        for (int i = 0; i < m_photonCount; ++i) {
            const Emitter * emitter = scene->getRandomEmitter(sampler->next1D());
            Ray3f ray;
            Color3f W = emitter->samplePhoton(ray, sampler->next2D(), sampler->next2D())*scene->getLights().size();
            tracePhoton(scene, ray, W, sampler);
        }

		/* Build the photon map */
        m_photonMap->build();
    }
        
    void tracePhoton(const Scene *scene, Ray3f &ray, Color3f &W, Sampler* sampler) {
        
        while (true) {
            
            // Check that photon bounces
            Intersection its;
            if (!scene->rayIntersect(ray, its)) {
                break;
            }
            
            // If diffuse surface, create photon and add it
            if (its.mesh->getBSDF()->isDiffuse()) {
                Photon photon(its.p, -ray.d, W);
                m_photonMap->push_back(photon);
            }
            
            // Russian roulette
            if (sampler->next1D() > std::min(W.maxCoeff(), 0.99f)) {
                break;
            }
            W /= std::min(W.maxCoeff(),0.999f);
            
            // Sample new direction
            BSDFQueryRecord bRec(its.toLocal(-ray.d));
            Color3f BSDF = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            ray = Ray3f(its.p, its.toWorld(bRec.wo));
            W *= BSDF;
        }
    }

    virtual Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const override {
    	
		/* How to find photons?
		 * std::vector<uint32_t> results;
		 * m_photonMap->search(Point3f(0, 0, 0), // lookup position
		 *                     m_photonRadius,   // search radius
		 *                     results);
		 *
		 * for (uint32_t i : results) {
		 *    const Photon &photon = (*m_photonMap)[i];
		 *    cout << "Found photon!" << endl;
		 *    cout << " Position  : " << photon.getPosition().toString() << endl;
		 *    cout << " Power     : " << photon.getPower().toString() << endl;
		 *    cout << " Direction : " << photon.getDirection().toString() << endl;
		 * }
		 */

		// put your code for path tracing with photon gathering here
        // Initial radiance and throughput
        Color3f Li(0.0f), t(1.f);
        Ray3f ray =_ray;
        
        while (true) {
            
            Intersection its;
            
            if(!scene->rayIntersect(ray, its)) {
                // No more intersection, return current Li
                return Li;
            }
            
            // First check if intersection is an emitter
            if (its.mesh->isEmitter()) {
                // Create a query to get the value
                EmitterQueryRecord lRec(ray.o, its.p, its.shFrame.n);
                Li += t*its.mesh->getEmitter()->eval(lRec);
            }
            
            // Contribution form photons
            if (its.mesh->getBSDF()->isDiffuse()) {
               return Li + t * photonDensityEstimation(its, ray);
            }
            
            // Russian roulette
            if (sampler->next1D() > std::min(t.maxCoeff(),0.999f)) {
                return Li;
            }
            t /= std::min(t.maxCoeff(),0.999f);
            
            // Use BSDF Sampling and shoot a ray in that direction
            BSDFQueryRecord bRec(its.shFrame.toLocal(-ray.d));
            bRec.uv = its.uv;
            Color3f BSDF = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            ray = Ray3f(its.p, its.toWorld(bRec.wo));
            
            // The sample function already returns the value divided by the pdf
            t *= BSDF;
            
        }
    }
        
    Color3f photonDensityEstimation(const Intersection &its, const Ray3f &ray) const {
        
        std::vector<Photon::IndexType> scopePhotons;
        m_photonMap->search(its.p, m_photonRadius, scopePhotons);
        
        Color3f Lr(0.0f);
        //BSDFQueryRecord bRec(its.shFrame.toLocal(-ray.d));
        
        for (int i = 0; i < scopePhotons.size(); ++i) {
            Photon photon = (*m_photonMap)[scopePhotons[i]];
            BSDFQueryRecord bRec(its.shFrame.toLocal(-ray.d), its.shFrame.toLocal(photon.getDirection()), ESolidAngle);
            Lr += its.mesh->getBSDF()->eval(bRec)*photon.getPower();
        }
            
        return Lr * INV_PI / (pow(m_photonRadius,2) * m_photonCount);
    }

    virtual std::string toString() const override {
        return tfm::format(
            "PhotonMapper[\n"
            "  photonCount = %i,\n"
            "  photonRadius = %f\n"
            "]",
            m_photonCount,
            m_photonRadius
        );
    }
private:
    int m_photonCount;
    float m_photonRadius;
    std::unique_ptr<PhotonMap> m_photonMap;
};

NORI_REGISTER_CLASS(PhotonMapper, "photonmapper");
NORI_NAMESPACE_END
