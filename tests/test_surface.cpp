#include "../tinynurbs/tinynurbs.h"
#include "glm/glm.hpp"
#include "glm/gtc/constants.hpp"
#include <cmath>

#include "catch.hpp"

using namespace std;

tinynurbs::RationalSurface3f getHemisphere() {
    tinynurbs::RationalSurface3f srf;
    srf.degree_u = 3;
    srf.degree_v = 3;
    srf.knots_u = {0, 0, 0, 0, 1, 1, 1, 1};
    srf.knots_v = {0, 0, 0, 0, 1, 1, 1, 1};
    // 4x4 grid (tinynurbs::array2) of control points and weights
    // https://www.geometrictools.com/Documentation/NURBSCircleSphere.pdf
    srf.control_points = {4, 4, 
                          {glm::vec3(0, 0, 1), glm::vec3(0, 0, 1), glm::vec3(0, 0, 1), glm::vec3(0, 0, 1),
                           glm::vec3(2, 0, 1), glm::vec3(2, 4, 1),  glm::vec3(-2, 4, 1),  glm::vec3(-2, 0, 1),
                           glm::vec3(2, 0, -1), glm::vec3(2, 4, -1), glm::vec3(-2, 4, -1), glm::vec3(-2, 0, -1),
                           glm::vec3(0, 0, -1), glm::vec3(0, 0, -1), glm::vec3(0, 0, -1), glm::vec3(0, 0, -1)
                          }
                         };
    srf.weights = {4, 4,
                   {1,       1.f/3.f, 1.f/3.f, 1,
                    1.f/3.f, 1.f/9.f, 1.f/9.f, 1.f/3.f,
                    1.f/3.f, 1.f/9.f, 1.f/9.f, 1.f/3.f,
                    1,       1.f/3.f, 1.f/3.f, 1
                   }
                 };
    return srf;
}

TEST_CASE("surfacePoint", "[surface, evaluate]")
{
    auto srf = getHemisphere();
    glm::vec3 pt1 = tinynurbs::surfacePoint(srf, 0.f, 0.f);
    glm::vec3 pt2 = tinynurbs::surfacePoint(srf, 1.f, 1.f);
    REQUIRE(pt1.x == Approx(0));
    REQUIRE(pt1.y == Approx(0));
    REQUIRE(pt1.z == Approx(1));
    REQUIRE(pt2.x == Approx(0));
    REQUIRE(pt2.y == Approx(0));
    REQUIRE(pt2.z == Approx(-1));
}

TEST_CASE("surfaceTangent", "[surface, evaluate]")
{
    auto srf = getHemisphere();
    glm::vec3 tgt_u, tgt_v;
    std::tie(tgt_u, tgt_v) = tinynurbs::surfaceTangent(srf, 0.f, 0.f);
    REQUIRE(glm::length(tgt_u) == Approx(1));
    REQUIRE(tgt_u.x == Approx(1));
    REQUIRE(tgt_u.y == Approx(0));
    REQUIRE(tgt_u.z == Approx(0));
    // tgt_v should be a zero vector since the pole of the hemisphere is pinched
    REQUIRE(glm::length(tgt_v) == Approx(0));
    REQUIRE(tgt_v.x == Approx(0));
    REQUIRE(tgt_v.y == Approx(0));
    REQUIRE(tgt_v.z == Approx(0));
}


TEST_CASE("surfaceNormal", "[surface, evaluate]")
{
    auto srf = getHemisphere();
    glm::vec3 n = tinynurbs::surfaceNormal(srf, 0.5f, 0.5f);
    REQUIRE(glm::length(n) == Approx(1));
    REQUIRE(n.x == Approx(0));
    REQUIRE(n.y == Approx(-1));
    REQUIRE(n.z == Approx(0));

    // Normal is zero vector since the pole is pinched
    n = tinynurbs::surfaceNormal(srf, 0.f, 0.f);
    REQUIRE(glm::length(n) == Approx(0));
    REQUIRE(n.x == Approx(0));
    REQUIRE(n.y == Approx(0));
    REQUIRE(n.z == Approx(0));
}

TEST_CASE("surfaceIsValid", "[surface, check]")
{
    auto srf = getHemisphere();
    bool is_valid = tinynurbs::surfaceIsValid(srf);
    REQUIRE(is_valid == true);

    srf = getHemisphere();
    srf.degree_u = 5;
    is_valid = tinynurbs::surfaceIsValid(srf);
    REQUIRE(is_valid == false);

    srf = getHemisphere();
    srf.degree_v = 4;
    is_valid = tinynurbs::surfaceIsValid(srf);
    REQUIRE(is_valid == false);
}