#include <nori/medium.h>

NORI_NAMESPACE_BEGIN

Color3f Medium::transmittance(const Point3f x, const Point3f y) const {
    float dist = (x - y).norm();
    return Color3f(exp(-m_sigmaT.x() * dist),
                   exp(-m_sigmaT.y() * dist),
                   exp(-m_sigmaT.z() * dist));
}

NORI_REGISTER_CLASS(Medium, "medium");
NORI_NAMESPACE_END
