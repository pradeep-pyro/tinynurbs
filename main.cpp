#include <iostream>
#define TINYNURBS_STATICLIB
#include "nurbstk/curve.h"
#include "nurbstk/evaluate.h"
#include "nurbstk/check.h"
#include "nurbstk/surface.h"
#include "nurbstk/io.h"

#define GLM_ENABLE_EXPERIMENTAL
#include "glm/glm.hpp"
#include "glm/gtx/string_cast.hpp"

using namespace std;

void testCurvePoint() {
    using namespace nurbstk;
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
    Curve<2, float> crv {degree, knots, controlPoints};
    glm::vec2 pt1, pt2, pt3;
    curvePoint(crv, 2.5f, pt1);
    curvePoint(crv, 0.f, pt2);
    curvePoint(crv, 5.f, pt3);
    cout << glm::to_string(pt1) << endl;
    cout << glm::to_string(pt2) << endl;
    cout << glm::to_string(pt3) << endl;
}

void testRationalCurvePoint() {
    using namespace nurbstk;
    unsigned int degree = 2;
    std::vector<float> knots {0, 0, 0, 1, 1, 1};
    std::vector<glm::vec2> controlPoints;
    controlPoints.push_back(glm::vec2(1, 0));
    controlPoints.push_back(glm::vec2(1, 1));
    controlPoints.push_back(glm::vec2(0, 1));
    RationalCurve2f crv {degree, knots, controlPoints, std::vector<float> {1, 1, 2}};
    glm::vec2 pt1, pt2, pt3;
    rationalCurvePoint(crv, 0.f, pt1);
    rationalCurvePoint(crv, 0.5f, pt2);
    rationalCurvePoint(crv, 1.f, pt3);
    cout << glm::to_string(pt1) << endl;
    cout << glm::to_string(pt2) << endl;
    cout << glm::to_string(pt3) << endl;
}

void testRationalSurfacePoint() {
    using namespace nurbstk;
    unsigned int degreeU = 1;
    unsigned int degreeV = 2;
    std::vector<float> knotsU {0, 0, 1, 1 }, knotsV {0, 0, 0, 1, 1, 1};
    nurbstk::array2<glm::vec3> cp;
    cp.resize(2, 3);
    cp(0, 0) = glm::vec3(1, 1, 0);
    cp(0, 1) = glm::vec3(1, 1, 1);
    cp(0, 2) = glm::vec3(1, 0, 1);
    cp(1, 0) = glm::vec3(-1, 1, 0);
    cp(1, 1) = glm::vec3(-1, 1, 1);
    cp(1, 2) = glm::vec3(-1, 0, 1);
    array2<float> w(2, 3);
    w(0, 0) = 1;
    w(0, 1) = 1;
    w(0, 2) = 2;
    w(1, 0) = 1;
    w(1, 1) = 1;
    w(1, 2) = 2;
    RationalSurface3f srf(degreeU, degreeV, knotsU, knotsV, cp, w);
    glm::vec3 pt1, pt2, pt3;
    rationalSurfacePoint(srf, 0.f, 0.f, pt1);
    rationalSurfacePoint(srf, 0.5f, 0.5f, pt2);
    rationalSurfacePoint(srf, 1.f, 1.f, pt3);
    cout << glm::to_string(pt1) << endl;
    cout << glm::to_string(pt2) << endl;
    cout << glm::to_string(pt3) << endl;
}

void testCurveDeriv() {
    using namespace nurbstk;
    unsigned int degree = 3;
    std::vector<float> knots {0, 0, 0, 0, 1, 1, 1, 1};
    std::vector<glm::vec2> controlPoints;
    controlPoints.push_back(glm::vec2(10, 0));
    controlPoints.push_back(glm::vec2(20, 10));
    controlPoints.push_back(glm::vec2(30, 20));
    controlPoints.push_back(glm::vec2(50, 50));
    nurbstk::Curve2f crv(degree, knots, controlPoints);
    std::vector<glm::vec2> ptder;
    curveDerivatives(crv, 2, 0.f, ptder);
    cout << "Derivative (order 0): " << ptder[0][0] << " " << ptder[0][1] << endl;
    cout << "Derivative (order 1): " << ptder[1][0] / ptder[1][1] << endl;
}

void testarray2() {
    std::vector<std::vector<float>> vec {{1, 1, 1}, {2, 2, 2}};
    nurbstk::array2<float> arr(2, 3);
    for (int i = 0; i < arr.rows(); ++i) {
        for (int j = 0; j < arr.cols(); ++j) {
            arr(i, j) = i == 0 ? 1 : 2;
        }
    }

    for (int i = 0; i < arr.rows(); ++i) {
        for (int j = 0; j < arr.cols(); ++j) {
            cout << vec[i][j] << " ";
        }
        cout << endl;
    }

    for (int i = 0; i < arr.rows(); ++i) {
        for (int j = 0; j < arr.cols(); ++j) {
            cout << arr(i, j) << " ";
        }
        cout << endl;
    }
}

void testCurveClosed() {
    unsigned int deg = 3;
    std::vector<float> knots {-0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75};
    std::vector<glm::vec2> cp;
    cp.push_back(glm::vec2(10, 0));
    cp.push_back(glm::vec2(20, 10));
    cp.push_back(glm::vec2(30, 20));
    cp.push_back(glm::vec2(50, 50));
    cp.push_back(glm::vec2(10, 0));
    cp.push_back(glm::vec2(20, 10));
    cp.push_back(glm::vec2(30, 20));

    nurbstk::Curve2f crv {deg, knots, cp};
    cout << nurbstk::isCurveClosed(crv) << endl;
}

int main() {
    testCurvePoint();
    testRationalCurvePoint();
    testCurveDeriv();
    testCurveClosed();
    testRationalSurfacePoint();
    return 0;
}
