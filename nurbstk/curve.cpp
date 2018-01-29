#include "curve.h"

namespace nurbstk {

// Explicit template class instantiations for
// 2D and 3D glm vector types
/*
template class NurbsCurve<2, float>;
template class NurbsCurve<2, double>;
template class NurbsCurve<3, float>;
template class NurbsCurve<3, double>;
*/
template class Curve<2, float>;
template class Curve<2, double>;
template class Curve<3, float>;
template class Curve<3, double>;

template class RationalCurve<2, float>;
template class RationalCurve<2, double>;
template class RationalCurve<3, float>;
template class RationalCurve<3, double>;

}