#include "../tinynurbs/tinynurbs.h"
#include "glm/glm.hpp"
#include "glm/gtc/constants.hpp"
#include <cmath>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace std;

// Unit circle
tinynurbs::RationalCurve3f getCircle() {
    tinynurbs::RationalCurve3f crv;
    crv.control_points = {glm::vec3(1, 0, 0),
                          glm::vec3(1, 1, 0),
                          glm::vec3(0, 1, 0),
                          glm::vec3(-1, 1, 0),
                          glm::vec3(-1, 0, 0),
                          glm::vec3(-1, -1, 0),
                          glm::vec3(0, -1, 0),
                          glm::vec3(1, -1, 0),
                          glm::vec3(1, 0, 0)
                         };
    const float sqrt2_over_2 = std::sqrt(2.f) / 2.f;
    crv.weights = {1, sqrt2_over_2, 1, sqrt2_over_2, 1,
                   sqrt2_over_2, 1, sqrt2_over_2, 1
                  };
    crv.knots = {0, 0, 0,
                 glm::half_pi<float>(), glm::half_pi<float>(),
                 glm::pi<float>(), glm::pi<float>(),
                 3 * glm::half_pi<float>(), 3 * glm::half_pi<float>(),
                 glm::two_pi<float>(), glm::two_pi<float>(), glm::two_pi<float>()
                };
    crv.degree = 2;
    return crv;
}

TEST_CASE("curvePoint", "[curve, evaluate]")
{
    auto crv = getCircle();
    glm::vec3 pt1;
    tinynurbs::curvePoint(crv, 0.f, pt1);
    REQUIRE(pt1.x == Approx(1));
    REQUIRE(pt1.y == Approx(0));
    REQUIRE(pt1.z == Approx(0));
    glm::vec3 pt2;
    curvePoint(crv, glm::pi<float>(), pt2);
    REQUIRE(pt2.x == Approx(-1));
    REQUIRE(pt2.y == Approx(0));
    REQUIRE(pt2.z == Approx(0));
}

TEST_CASE("curveTangent", "[curve, evaluate]")
{
    auto crv = getCircle();
    glm::vec3 tgt1;
    tinynurbs::curveTangent(crv, 0.f, tgt1);
    REQUIRE(tgt1.x == Approx(0));
    REQUIRE(tgt1.y == Approx(1));
    REQUIRE(tgt1.z == Approx(0));
    glm::vec3 tgt2;
    tinynurbs::curveTangent(crv, glm::pi<float>(), tgt2);
    REQUIRE(tgt2.x == Approx(0));
    REQUIRE(tgt2.y == Approx(-1));
    REQUIRE(tgt2.z == Approx(0));
}

TEST_CASE("curveIsValid", "[curve, check]")
{
    auto crv = getCircle();
    bool is_valid = tinynurbs::curveIsValid(crv);
    REQUIRE(is_valid == true);

    crv = getCircle();
    crv.degree = 4;
    is_valid = tinynurbs::curveIsValid(crv);
    REQUIRE(is_valid == false);
}

TEST_CASE("curveInsertKnot", "[curve, modify]")
{
    auto crv = getCircle();
    glm::vec3 pt, new_pt;
    tinynurbs::curvePoint(crv, 0.25f, pt);
    size_t n_knots_prev = crv.knots.size();
    size_t n_control_points_prev = crv.control_points.size();

    tinynurbs::curveKnotInsert(crv, 0.25f, 1);
    tinynurbs::curvePoint(crv, 0.25f, new_pt);
    size_t n_knots_curr = crv.knots.size();
    size_t n_control_points_curr = crv.control_points.size();

    REQUIRE((n_knots_prev + 1) == n_knots_curr);
    REQUIRE((n_control_points_prev + 1) == n_control_points_curr);
    REQUIRE(pt.x == Approx(new_pt.x));
    REQUIRE(pt.y == Approx(new_pt.y));
    REQUIRE(pt.z == Approx(new_pt.z));
}

TEST_CASE("curveSplit", "[curve, modify]")
{
    auto crv = getCircle();
    float u = glm::pi<float>();
    tinynurbs::RationalCurve3f left, right;
    tinynurbs::curveSplit(crv, u, left, right);

    bool is_valid_l = tinynurbs::curveIsValid(left);
    bool is_valid_r = tinynurbs::curveIsValid(right);

    REQUIRE(is_valid_l == true);
    REQUIRE(is_valid_r == true);

    REQUIRE(left.degree == crv.degree);
    REQUIRE(right.degree == crv.degree);

    for (int i = 0; i < left.degree + 1; ++i) {
        int d = left.knots.size() - (left.degree + 1);
        REQUIRE(left.knots[d+i] == Approx(u));
    }

    for (int i = 0; i < right.degree + 1; ++i) {
        REQUIRE(right.knots[i] == Approx(u));
    }

    glm::vec3 pt1, pt2;
    tinynurbs::curvePoint(crv, left.knots[left.knots.size() - 1], pt1);
    tinynurbs::curvePoint(crv, right.knots[0], pt2);
    REQUIRE(pt1.x == Approx(pt2.x));
    REQUIRE(pt1.y == Approx(pt2.y));
    REQUIRE(pt1.z == Approx(pt2.z));
}

/*
TEST_CASE("curveReadOBJ and curveSaveOBJ", "[curve, obj]")
{
    tinynurbs::RationalCurve3f crv;
    tinynurbs::curveReadOBJ("curve.obj", crv);
    REQUIRE(crv.degree == 2);
    REQUIRE(crv.knots.size() == 9);
    REQUIRE(crv.control_points.size() == 6);
    REQUIRE(crv.control_points[0].x == Approx(-2.3));
    REQUIRE(crv.control_points[3].y == Approx(-1.49));
    REQUIRE(crv.control_points[5].z == Approx(0));
    for (int i = 0; i < crv.control_points.size(); ++i) {
        REQUIRE(crv.weights[i] == Approx(1));
    }
    tinynurbs::curveSaveOBJ("curve_out.obj", crv);
}
*/