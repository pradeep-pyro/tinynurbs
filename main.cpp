#include <iostream>
#define TINYNURBS_STATICLIB
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
    nurbstk::array2<glm::vec3> cp;
    cp.resize(2, 3);
    cp(0, 0) = glm::vec3(1, 1, 0);
    cp(0, 1) = glm::vec3(1, 1, 1);
    cp(0, 2) = glm::vec3(1, 0, 1);
    cp(1, 0) = glm::vec3(-1, 1, 0);
    cp(1, 1) = glm::vec3(-1, 1, 1);
    cp(1, 2) = glm::vec3(-1, 0, 1);
    nurbstk::array2<float> w(2, 3, {1, 1, 2, 1, 1, 2});
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

    nurbstk::Curve2f crv(deg, knots, cp);
    cout << crv.isClosed() << endl;
}

int main() {
    testCurvePoint();
    testRationalCurvePoint();
    testCurveDeriv();
    testCurveClosed();
    return 0;
}
