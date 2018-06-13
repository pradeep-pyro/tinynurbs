#include "curve.h"

namespace nurbstk {

bool isValidRelation(unsigned int degree, size_t num_knots, size_t num_ctrl_pts) {
    return (num_knots - degree - 1) == num_ctrl_pts;
}

// Explicit template class instantiations for
// 2D and 3D glm vector types
template class Curve<2, float>;
template class Curve<2, double>;
template class Curve<3, float>;
template class Curve<3, double>;

template class RationalCurve<2, float>;
template class RationalCurve<2, double>;
template class RationalCurve<3, float>;
template class RationalCurve<3, double>;
}
