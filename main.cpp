#include "nurbstk/surface.h"
#include "nurbstk/io.h"

#include "glm/glm.hpp"

template <typename T>
struct SurfaceData {
    typedef glm::vec<4, T> vec;
    unsigned int deg_u, deg_v;
    std::vector<T> knots_u, knots_v;
    std::vector<vec> ctrl_pts;
    bool rational;
};

int main() {
    SurfaceData<float> data; // 3d bspline surface
    // SurfaceData<float, 4> data; // 4d rat bspline surface
    readOBJ("test.obj", data.deg_u, data.deg_v, data.knots_u, data.knots_v, data.ctrl_pts, data.rational);

    nurbstk::RationalSurface3f surf(data.deg_u, data.deg_v, data.knots_u, data.knots_v, data.ctrl_pts);
    surf.controlPoint(0, 0);

    return 0;
}
