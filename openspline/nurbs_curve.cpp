#include "nurbs_curve.h"

namespace st {

// Explicit template class instantiations for
// 2D and 3D glm vector types
template class NurbsCurve<2, float>;
template class NurbsCurve<2, double>;
template class NurbsCurve<3, float>;
template class NurbsCurve<3, double>;

}