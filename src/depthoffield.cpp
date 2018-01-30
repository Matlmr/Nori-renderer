/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob, Romain Prévost

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

#include <nori/camera.h>
#include <nori/rfilter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

/**
 * \brief Perspective camera with depth of field
 *
 * This class implements a perspective camera model with depth of field.
 */
class DepthOfFieldCamera : public Camera {
public:
    DepthOfFieldCamera(const PropertyList &propList) {
        /* Width and height in pixels. Default: 720p */
        m_outputSize.x() = propList.getInteger("width", 1280);
        m_outputSize.y() = propList.getInteger("height", 720);
        m_invOutputSize = m_outputSize.cast<float>().cwiseInverse();

        /* Specifies an optional camera-to-world transformation. Default: none */
        m_cameraToWorld = propList.getTransform("toWorld", Transform());

        /* Horizontal field of view in degrees */
        m_fov = propList.getFloat("fov", 30.0f);

        /* Near and far clipping planes in world-space units */
        m_nearClip = propList.getFloat("nearClip", 1e-4f);
        m_farClip = propList.getFloat("farClip", 1e4f);
        
        /* Lens radius and focal distance needed for Depth of Field */
        m_lensRadius = propList.getFloat("lensradius", -1);
        m_focalDistance = propList.getFloat("focaldistance");

        m_rfilter = NULL;
    }

    virtual void activate() override {
        float aspect = m_outputSize.x() / (float) m_outputSize.y();

        /* Project vectors in camera space onto a plane at z=1:
         *
         *  xProj = cot * x / z
         *  yProj = cot * y / z
         *  zProj = (far * (z - near)) / (z * (far-near))
         *  The cotangent factor ensures that the field of view is 
         *  mapped to the interval [-1, 1].
         */
        float recip = 1.0f / (m_farClip - m_nearClip),
              cot = 1.0f / std::tan(degToRad(m_fov / 2.0f));

        Eigen::Matrix4f perspective;
        perspective <<
            cot, 0,   0,   0,
            0, cot,   0,   0,
            0,   0,   m_farClip * recip, -m_nearClip * m_farClip * recip,
            0,   0,   1,   0;

        /**
         * Translation and scaling to shift the clip coordinates into the
         * range from zero to one. Also takes the aspect ratio into account.
         */
        m_sampleToCamera = Transform( 
            Eigen::DiagonalMatrix<float, 3>(Vector3f(0.5f, -0.5f * aspect, 1.0f)) *
            Eigen::Translation<float, 3>(1.0f, -1.0f/aspect, 0.0f) * perspective).inverse();

        /* If no reconstruction filter was assigned, instantiate a Gaussian filter */
        if (!m_rfilter) {
            m_rfilter = static_cast<ReconstructionFilter *>(
                    NoriObjectFactory::createInstance("gaussian", PropertyList()));
            m_rfilter->activate();
        }
    }

    Color3f sampleRay(Ray3f &ray,
            const Point2f &samplePosition,
            const Point2f &apertureSample) const {
        /* Compute the corresponding position on the 
           near plane (in local camera space) */
        Point3f nearP = m_sampleToCamera * Point3f(
            samplePosition.x() * m_invOutputSize.x(),
            samplePosition.y() * m_invOutputSize.y(), 0.0f);

        /* Turn into a normalized ray direction, and
           adjust the ray interval accordingly */
        Vector3f d = nearP.normalized();
        float invZ = 1.0f / d.z();
        
        if (m_lensRadius > 0.f) {
            // Sample point on lens
            Point2f lens = Warp::squareToUniformDisk(apertureSample);
            //Point2f lens = ConcentricSampleDisk(apertureSample);
            float lensU = lens.x() * m_lensRadius;
            float lensV = lens.y() * m_lensRadius;
            // Compute point on plane of focus
            float ft = m_focalDistance / d.z();
            Point3f Pfocus = d*ft;
            // Update ray for effect of lens
            ray.o = m_cameraToWorld * Point3f(lensU, lensV, 0);
            ray.d = m_cameraToWorld * Vector3f((Pfocus - Point3f(lensU, lensV,0)).normalized());
        } else {
            // Like for perspective camera
            ray.o = m_cameraToWorld * Point3f(0, 0, 0);
            ray.d = m_cameraToWorld * d;
        }
        ray.mint = m_nearClip * invZ;
        ray.maxt = m_farClip * invZ;
        ray.update();
        
        return Color3f(1.0f);
    }
    
    Point2f ConcentricSampleDisk(const Point2f &sample) const {
        float r, theta;
        float sx = 2 * sample.x() - 1;
        float sy = 2 * sample.y() - 1;
        if (sx >= -sy) {
            if (sx > sy) {
                // First region of disk
                r = sx;
                if (sy > 0.0f) theta = sy/r;
                else theta = 8.0f + sy/r;
            } else {
                // Second region
                if (sx > 0.0f) r = sx;
                else r = - sx;
                theta = sy/r;
            }
        } else {
            if (sx <= sy) {
                // Third region
                r = - sx;
                if (sy > 0.0f) theta = sy/r;
                else theta = 8.0f + sy/r;
            } else {
                // Fourth region
                if (sx > 0.0f) r = sx;
                else r = - sx;
                theta = 8.0f + sy/r;
            }
        }
        theta *= M_PI / 4.0f;
        return Point2f(r * cos(theta), r * sin(theta));
    }
    

    virtual void addChild(NoriObject *obj) override {
        switch (obj->getClassType()) {
            case EReconstructionFilter:
                if (m_rfilter)
                    throw NoriException("Camera: tried to register multiple reconstruction filters!");
                m_rfilter = static_cast<ReconstructionFilter *>(obj);
                break;

            default:
                throw NoriException("Camera::addChild(<%s>) is not supported!",
                    classTypeName(obj->getClassType()));
        }
    }

    /// Return a human-readable summary
    virtual std::string toString() const override {
        return tfm::format(
            "PerspectiveCamera[\n"
            "  cameraToWorld = %s,\n"
            "  outputSize = %s,\n"
            "  fov = %f,\n"
            "  clip = [%f, %f],\n"
            "  lensRadius = %f,\n"
            "  focalDistance = %f,\n"
            "  rfilter = %s\n"
            "]",
            indent(m_cameraToWorld.toString(), 18),
            m_outputSize.toString(),
            m_fov,
            m_nearClip,
            m_farClip,
            m_lensRadius,
            m_focalDistance,
            indent(m_rfilter->toString())
        );
    }
private:
    Vector2f m_invOutputSize;
    Transform m_sampleToCamera;
    Transform m_cameraToWorld;
    float m_fov;
    float m_nearClip;
    float m_farClip;
    float m_lensRadius;
    float m_focalDistance;
};

NORI_REGISTER_CLASS(DepthOfFieldCamera, "depthoffield");
NORI_NAMESPACE_END
