#include <iostream>
#include "tinynurbs/geometry/curve.h"
#include "tinynurbs/core/evaluate.h"
#include "tinynurbs/core/check.h"
#include "tinynurbs/geometry/surface.h"
#include "tinynurbs/io/obj.h"
#include "tinynurbs/modify/refine.h"

#define GLM_ENABLE_EXPERIMENTAL
#include "glm/glm.hpp"
#include "glm/gtx/string_cast.hpp"
#include <cmath>

using namespace std;

void testCurvePoint() {
    using namespace tinynurbs;
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
    using namespace tinynurbs;
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
    using namespace tinynurbs;
    unsigned int degreeU = 1;
    unsigned int degreeV = 2;
    std::vector<float> knotsU {0, 0, 1, 1 }, knotsV {0, 0, 0, 1, 1, 1};
    array2<glm::vec3> cp;
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
    using namespace tinynurbs;
    unsigned int degree = 3;
    std::vector<float> knots {0, 0, 0, 0, 1, 1, 1, 1};
    std::vector<glm::vec2> controlPoints;
    controlPoints.push_back(glm::vec2(10, 0));
    controlPoints.push_back(glm::vec2(20, 10));
    controlPoints.push_back(glm::vec2(30, 20));
    controlPoints.push_back(glm::vec2(50, 50));
    Curve2f crv(degree, knots, controlPoints);
    std::vector<glm::vec2> ptder;
    curveDerivatives(crv, 2, 0.f, ptder);
    cout << "Derivative (order 0): " << ptder[0][0] << " " << ptder[0][1] << endl;
    cout << "Derivative (order 1): " << ptder[1][0] / ptder[1][1] << endl;
}

void testCurveClosed() {
    using namespace tinynurbs;
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

    Curve2f crv {deg, knots, cp};
    cout << isCurveClosed(crv) << endl;
}

void testIO() {
    unsigned int degu, degv;
    std::vector<float> knotsu, knotsv;
    tinynurbs::array2<glm::vec3> cp;
    tinynurbs::array2<float> w;
    bool rat;
    tinynurbs::readOBJ("/home/pradeep/Downloads/car50_100.obj", degu, degv, knotsu, knotsv, cp, w, rat);
    tinynurbs::saveOBJ("/home/pradeep/Downloads/car50_1002.obj", degu, degv, knotsu, knotsv, cp, w, rat);
}

void testKnotInsert() {
    using namespace tinynurbs;
    Curve2f crv;
    crv.degree = 3;
    crv.control_points = {{5.0, 5.0}, {10.0, 10.0}, {20.0, 15.0}, {35.0, 15.0}, {45.0, 10.0}, {50.0, 5.0}};
    crv.knots = {0.0, 0.0, 0.0, 0.0, 0.33, 0.66, 1.0, 1.0, 1.0, 1.0};

    // Set evaluation parameter
    float u = 0.3f;
    glm::vec2 pt;
    curvePoint(crv, u, pt);
    cout << glm::to_string(pt) << endl;
    curveInsertKnot(crv.degree, crv.knots, crv.control_points, u);
    curvePoint(crv, u, pt);
    cout << glm::to_string(pt) << endl;
    // Evaluation result
    glm::vec2 res(18.617, 13.377);
}

int main() {
    /*testCurvePoint();
    testRationalCurvePoint();
    testCurveDeriv();
    testCurveClosed();
    testRationalSurfacePoint();
    testIO();*/
    testKnotInsert();
    return 0;
}
