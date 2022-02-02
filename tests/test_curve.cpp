#include <tinynurbs/tinynurbs.h>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <cmath>

#include "catch.hpp"

using namespace std;

tinynurbs::Curve3f getNonrationalBezierCurve() {
    tinynurbs::Curve3f crv;
    crv.control_points = {glm::vec3(-1, 0, 0),
                          glm::vec3(0, 1, 0),
                          glm::vec3(1, 0, 0)
                         };
    crv.knots = {0, 0, 0, 1, 1, 1};
    crv.degree = 2;
    return crv;
}

TEST_CASE("curvePoint (non-rational)", "[curve, non-rational, evaluate]")
{
    auto crv = getNonrationalBezierCurve();
    glm::vec2 pt1 = tinynurbs::curvePoint(crv, 0.f);
    REQUIRE(pt1.x == Approx(-1));
    REQUIRE(pt1.y == Approx(0));
    glm::vec2 pt2 = tinynurbs::curvePoint(crv, 1.f);
    REQUIRE(pt2.x == Approx(1));
    REQUIRE(pt2.y == Approx(0));
}

TEST_CASE("curveTangent (glm::vec2)", "[curve, non-rational, evaluate]")
{
    auto crv = getNonrationalBezierCurve();
    glm::vec2 tgt1 = tinynurbs::curveTangent(crv, 0.5f);
    REQUIRE(tgt1.x == Approx(1));
    REQUIRE(tgt1.y == Approx(0));
}

TEST_CASE("curveIsValid (non-rational)", "[curve, non-rational, check]")
{
    auto crv = getNonrationalBezierCurve();
    bool is_valid = tinynurbs::curveIsValid(crv);
    REQUIRE(is_valid == true);

    crv = getNonrationalBezierCurve();
    crv.degree = 4;
    is_valid = tinynurbs::curveIsValid(crv);
    REQUIRE(is_valid == false);
}

TEST_CASE("curveKnotMultiplicity (non-rational)", "[knots, check]")
{
    {
        const std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 1.0};
        unsigned int knotMult0 = tinynurbs::knotMultiplicity(knots, 0.0);
        REQUIRE(knotMult0 == 3);
        unsigned int knotMult1 = tinynurbs::knotMultiplicity(knots, 1.0);
        REQUIRE(knotMult1 == 2);
    }
    {
        const std::vector<double> knots{0.0, 1.0};
        unsigned int knotMult0 = tinynurbs::knotMultiplicity(knots, 0.0);
        REQUIRE(knotMult0 == 1);
        unsigned int knotMult1 = tinynurbs::knotMultiplicity(knots, 1.0);
        REQUIRE(knotMult1 == 1);
    }
    {
        const std::vector<double> knots{0.0, 0.0001, 0.0002, 1.0, 1.0001};
        unsigned int knotMult0 = tinynurbs::knotMultiplicity(knots, 0.0);
        REQUIRE(knotMult0 == 1);
        unsigned int knotMult1 = tinynurbs::knotMultiplicity(knots, 1.0);
        REQUIRE(knotMult1 == 1);
    }
}

TEST_CASE("curveInsertKnot (non-rational)", "[curve, non-rational, modify]")
{
    auto crv = getNonrationalBezierCurve();
    glm::vec2 pt = tinynurbs::curvePoint(crv, 0.25f);
    
    size_t n_knots_prev = crv.knots.size();
    size_t n_control_points_prev = crv.control_points.size();

    auto new_crv = tinynurbs::curveKnotInsert(crv, 0.25f, 1);
    glm::vec2 new_pt = tinynurbs::curvePoint(new_crv, 0.25f);
    size_t n_knots_curr = new_crv.knots.size();
    size_t n_control_points_curr = new_crv.control_points.size();

    REQUIRE((n_knots_prev + 1) == n_knots_curr);
    REQUIRE((n_control_points_prev + 1) == n_control_points_curr);
    REQUIRE(pt.x == Approx(new_pt.x));
    REQUIRE(pt.y == Approx(new_pt.y));
}

TEST_CASE("curveSplit (non-rational)", "[curve, non-rational, modify]")
{
    auto crv = getNonrationalBezierCurve();
    float u = 0.5f;
    tinynurbs::Curve3f left, right;
    std::tie(left, right) = tinynurbs::curveSplit(crv, u);

    bool is_valid_l = tinynurbs::curveIsValid(left);
    bool is_valid_r = tinynurbs::curveIsValid(right);

    REQUIRE(is_valid_l == true);
    REQUIRE(is_valid_r == true);

    REQUIRE(left.degree == crv.degree);
    REQUIRE(right.degree == crv.degree);

    for (unsigned int i = 0; i < left.degree + 1; ++i) {
        int d = static_cast<int>(left.knots.size()) - (left.degree + 1);
        REQUIRE(left.knots[d+i] == Approx(u));
    }

    for (unsigned int i = 0; i < right.degree + 1; ++i) {
        REQUIRE(right.knots[i] == Approx(u));
    }

    glm::vec2 pt1 = tinynurbs::curvePoint(crv, left.knots[left.knots.size() - 1]);
    glm::vec2 pt2 = tinynurbs::curvePoint(crv, right.knots[0]);
    REQUIRE(pt1.x == Approx(pt2.x));
    REQUIRE(pt1.y == Approx(pt2.y));
}

TEST_CASE("curveReadOBJ and curveSaveOBJ (non-rational)", "[curve, non-rational, obj]")
{
    auto crv = getNonrationalBezierCurve();
    tinynurbs::curveSaveOBJ("curve.obj", crv);
    auto read_crv = tinynurbs::Curve3f(tinynurbs::curveReadOBJ<float>("curve.obj"));
    REQUIRE(crv.degree == read_crv.degree);
    REQUIRE(crv.knots.size() == read_crv.knots.size());
    for (int i = 0; i < crv.knots.size(); ++i) {
        REQUIRE(crv.knots[i] == Approx(read_crv.knots[i]));
    }
    REQUIRE(crv.control_points.size() == read_crv.control_points.size());
    for (int i = 0; i < crv.control_points.size(); ++i) {
        REQUIRE(crv.control_points[i].x == Approx(read_crv.control_points[i].x));
        REQUIRE(crv.control_points[i].y == Approx(read_crv.control_points[i].y));
    }
}