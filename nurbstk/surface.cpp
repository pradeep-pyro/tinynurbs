#include "surface.h"

#ifdef TINYNURBS_STATICLIB
namespace tinynurbs {

template class Surface<3, float>;
template class Surface<3, double>;
template class RationalSurface<3, float>;
template class RationalSurface<3, double>;

} // namespace nurbstk
#endif