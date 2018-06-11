#include <iostream>
#include "nurbstk/curve.h"
#include "nurbstk/surface.h"
#include "nurbstk/io.h"

#define GLM_ENABLE_EXPERIMENTAL

#include "glm/glm.hpp"
#include "glm/gtx/string_cast.hpp"

using namespace std;

void testCurvePoint() {
    unsigned int degree = 2;
    std::vector<float> knots = {0, 0, 0, 1, 2, 3, 4, 5, 5, 5};
    std::vector<glm::vec2> controlPoints;
    controlPoints.push_back(glm::vec2(10, 0));
    controlPoints.push_back(glm::vec2(20, 10));
    controlPoints.push_back(glm::vec2(30, 20));
    controlPoints.push_back(glm::vec2(40, 30));
    controlPoints.push_back(glm::vec2(50, 40));
    controlPoints.push_back(glm::vec2(60, 30));
    controlPoints.push_back(glm::vec2(70, 80));
    nurbstk::Curve2f crv(degree, knots, controlPoints);
    cout << glm::to_string(crv.point(2.5)) << endl;
    cout << glm::to_string(crv.point(0.)) << endl;
    cout << glm::to_string(crv.point(5.)) << endl;
}

void testRationalCurvePoint() {
    unsigned int degree = 2;
    std::vector<float> knots {0, 0, 0, 1, 1, 1};
    std::vector<glm::vec2> controlPoints;
    controlPoints.push_back(glm::vec2(1, 0));
    controlPoints.push_back(glm::vec2(1, 1));
    controlPoints.push_back(glm::vec2(0, 1));
    nurbstk::RationalCurve2f crv(degree, knots, controlPoints, std::vector<float> {1, 1, 2});
    cout << glm::to_string(crv.point(0)) << endl;
    cout << glm::to_string(crv.point(0.5)) << endl;
    cout << glm::to_string(crv.point(1)) << endl;
}

void testRationalSurfacePoint() {
    unsigned int degreeU = 1;
    unsigned int degreeV = 2;
    std::vector<float> knotsU {0, 0, 1, 1 }, knotsV {0, 0, 0, 1, 1, 1};
    std::vector<std::vector<glm::vec3>> cp;
    cp.resize(2);
    cp[0].push_back(glm::vec3(1, 1, 0));
    cp[0].push_back(glm::vec3(1, 1, 1));
    cp[0].push_back(glm::vec3(1, 0, 1));
    cp[1].push_back(glm::vec3(-1, 1, 0));
    cp[1].push_back(glm::vec3(-1, 1, 1));
    cp[1].push_back(glm::vec3(-1, 0, 1));
    std::vector<std::vector<float>> w {{1, 1, 2}, {1, 1, 2}};
    nurbstk::RationalSurface3f srf(degreeU, degreeV, knotsU, knotsV, cp, w);
    cout << glm::to_string(srf.point(0, 0)) << endl;
    cout << glm::to_string(srf.point(0.5, 0.5)) << endl;
    cout << glm::to_string(srf.point(1, 1)) << endl;
}

void testCurveDeriv() {
    unsigned int degree = 3;
    std::vector<float> knots {0, 0, 0, 0, 1, 1, 1, 1};
    std::vector<glm::vec2> controlPoints;
    controlPoints.push_back(glm::vec2(10, 0));
    controlPoints.push_back(glm::vec2(20, 10));
    controlPoints.push_back(glm::vec2(30, 20));
    controlPoints.push_back(glm::vec2(50, 50));
    nurbstk::Curve2f crv(degree, knots, controlPoints);
    std::vector<glm::vec2> ptder;
    crv.derivatives(0, 2, ptder);
    cout << "Derivative (order 0): " << glm::to_string(ptder[0]) << endl;
    cout << "Derivative (order 1): " << glm::to_string(ptder[1]) << endl;
    cout << "Derivative (order 2): " << glm::to_string(ptder[2]) << endl;
}

int main() {
    nurbstk::Surface3d srf = nurbstk::Surface3d::fromOBJ("../car50_100.obj");
    cout << glm::to_string(srf.point(0.5, 0.5)) << endl;

    unsigned int degreeU;
    unsigned int degreeV;
    std::vector<float> knotsU, knotsV;
    nurbstk::util::array2<glm::vec3> cp;
    nurbstk::util::array2<float> w;
    bool rational;
    nurbstk::readOBJ<float>("../car50_100.obj", degreeU, degreeV, knotsU, knotsV, cp, w, rational);
    for (int i = 0; i < cp.size(); ++i) {
        cout << glm::to_string(cp[i]) << endl;
    }
    return 0;
}
