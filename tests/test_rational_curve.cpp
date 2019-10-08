#include <tinynurbs/tinynurbs.h>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <cmath>

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

TEST_CASE("curvePoint (rational)", "[curve, rational, evaluate]")
{
    auto crv = getCircle();
    glm::vec3 pt1 = tinynurbs::curvePoint(crv, 0.f);
    REQUIRE(pt1.x == Approx(1));
    REQUIRE(pt1.y == Approx(0));
    REQUIRE(pt1.z == Approx(0));
    glm::vec3 pt2 = curvePoint(crv, glm::pi<float>());
    REQUIRE(pt2.x == Approx(-1));
    REQUIRE(pt2.y == Approx(0));
    REQUIRE(pt2.z == Approx(0));
}

TEST_CASE("curveTangent (rational)", "[curve, rational, evaluate]")
{
    auto crv = getCircle();
    glm::vec3 tgt1 = tinynurbs::curveTangent(crv, 0.f);
    REQUIRE(tgt1.x == Approx(0));
    REQUIRE(tgt1.y == Approx(1));
    REQUIRE(tgt1.z == Approx(0));
    glm::vec3 tgt2 = tinynurbs::curveTangent(crv, glm::pi<float>());
    REQUIRE(tgt2.x == Approx(0));
    REQUIRE(tgt2.y == Approx(-1));
    REQUIRE(tgt2.z == Approx(0));
}

TEST_CASE("curveIsValid (rational)", "[curve, rational, check]")
{
    auto crv = getCircle();
    bool is_valid = tinynurbs::curveIsValid(crv);
    REQUIRE(is_valid == true);

    crv = getCircle();
    crv.degree = 4;
    is_valid = tinynurbs::curveIsValid(crv);
    REQUIRE(is_valid == false);
}

TEST_CASE("curveInsertKnot (rational)", "[curve, rational, modify]")
{
    auto crv = getCircle();
    glm::vec3 pt = tinynurbs::curvePoint(crv, 0.25f);
    
    size_t n_knots_prev = crv.knots.size();
    size_t n_control_points_prev = crv.control_points.size();

    auto new_crv = tinynurbs::curveKnotInsert(crv, 0.25f, 1);
    glm::vec3 new_pt = tinynurbs::curvePoint(new_crv, 0.25f);
    size_t n_knots_curr = new_crv.knots.size();
    size_t n_control_points_curr = new_crv.control_points.size();

    REQUIRE((n_knots_prev + 1) == n_knots_curr);
    REQUIRE((n_control_points_prev + 1) == n_control_points_curr);
    REQUIRE(pt.x == Approx(new_pt.x));
    REQUIRE(pt.y == Approx(new_pt.y));
    REQUIRE(pt.z == Approx(new_pt.z));
}

TEST_CASE("curveSplit (rational)", "[curve, rational, modify]")
{
    auto crv = getCircle();
    float u = glm::pi<float>();
    tinynurbs::RationalCurve3f left, right;
    std::tie(left, right) = tinynurbs::curveSplit(crv, u);

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

    glm::vec3 pt1 = tinynurbs::curvePoint(crv, left.knots[left.knots.size() - 1]);
    glm::vec3 pt2 = tinynurbs::curvePoint(crv, right.knots[0]);
    REQUIRE(pt1.x == Approx(pt2.x));
    REQUIRE(pt1.y == Approx(pt2.y));
    REQUIRE(pt1.z == Approx(pt2.z));
}

TEST_CASE("curveReadOBJ and curveSaveOBJ (rational)", "[curve, rational, obj]")
{
    auto crv = getCircle();
    tinynurbs::curveSaveOBJ("curve.obj", crv);
    auto read_crv = tinynurbs::curveReadOBJ<float>("curve.obj");
    REQUIRE(crv.degree == read_crv.degree);
    REQUIRE(crv.knots.size() == read_crv.knots.size());
    for (int i = 0; i < crv.knots.size(); ++i) {
        REQUIRE(crv.knots[i] == Approx(read_crv.knots[i]));
    }
    REQUIRE(crv.control_points.size() == read_crv.control_points.size());
    for (int i = 0; i < crv.control_points.size(); ++i) {
        REQUIRE(crv.control_points[i].x == Approx(read_crv.control_points[i].x));
        REQUIRE(crv.control_points[i].y == Approx(read_crv.control_points[i].y));
        REQUIRE(crv.control_points[i].z == Approx(read_crv.control_points[i].z));
    }
    REQUIRE(crv.weights.size() == read_crv.weights.size());
    for (int i = 0; i < crv.weights.size(); ++i) {
        REQUIRE(crv.weights[i] == Approx(read_crv.weights[i]));
    }
}