#include <iostream>
#include "nurbstk/curve.h"

#define GLM_ENABLE_EXPERIMENTAL

#include "glm/glm.hpp"
#include "glm/gtx/string_cast.hpp"

using namespace std;

template <typename T>
struct SurfaceData {
    typedef glm::vec<4, T> vec;
    unsigned int deg_u, deg_v;
    std::vector<T> knots_u, knots_v;
    std::vector<vec> ctrl_pts;
    bool rational;
};

int main() {
    // SurfaceData<float> data; // 3d bspline surface
    // SurfaceData<float, 4> data; // 4d rat bspline surface
    // readOBJ("test.obj", data.deg_u, data.deg_v, data.knots_u, data.knots_v, data.ctrl_pts, data.rational);

    nurbstk::Curve2f curve(3, std::vector<float> {0.f, 0.f, 0.f, 0.f, 1.f, 1.f, 1.f, 1.f}, std::vector<glm::vec2> {glm::vec2(0, 0), glm::vec2(1, 0), glm::vec2(1, 1), glm::vec2(2, 0)});
    cout << glm::to_string(curve.point(0.5)) << endl;
    return 0;
}
