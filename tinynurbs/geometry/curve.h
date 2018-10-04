/**
@file
@brief The Curve class represents a non-uniform polynomial B-spline curve, while the RationalCurve class
represents a non-uniform rational B-spline (NURBS) curve.

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE file.
*/

#ifndef TINYNURBS_CURVE_H
#define TINYNURBS_CURVE_H

#include <vector>
#include <exception>
#include <stdexcept>
#include "glm/glm.hpp"

namespace tinynurbs {

// Forward declaration
template <int dim, typename T>
struct RationalCurve;

/**
Struct for holding a polynomial B-spline curve
@tparam dim Dimension of the curve (2 or 3)
@tparam T Data type of control points and knots (float or double)
*/
template <int dim, typename T>
struct Curve {
    unsigned int degree;
    std::vector<T> knots;
    std::vector<glm::vec<dim, T>> control_points;

    Curve() = default;
    Curve(const RationalCurve<dim, T> &crv) : Curve(crv.degree, crv.knots, crv.control_points) {
    }
    Curve(unsigned int degree, const std::vector<T> &knots,
          const std::vector<glm::vec<dim, T>> &control_points)
        : degree(degree), knots(knots), control_points(control_points) {
    }
};

/**
Struct for holding a rational B-spline curve
@tparam dim Dimension of the curve (2 or 3)
@tparam T Data type of control points and knots (float or double)
*/
template <int dim, typename T>
struct RationalCurve {
    unsigned int degree;
    std::vector<T> knots;
    std::vector<glm::vec<dim, T>> control_points;
    std::vector<T> weights;

    RationalCurve() = default;
    RationalCurve(const Curve<dim, T> &crv)
        : RationalCurve(crv, std::vector<T>(crv.control_points.size(), 1.0)) {
    }
    RationalCurve(const Curve<dim, T> &crv, const std::vector<T> &weights)
        : RationalCurve(crv.degree, crv.knots, crv.control_points, weights) {
    }
    RationalCurve(unsigned int degree, const std::vector<T> &knots,
                  const std::vector<glm::vec<dim, T>> &control_points,
                  const std::vector<T> weights)
        : degree(degree), knots(knots), control_points(control_points),
          weights(weights) {
    }
};

// Typedefs for ease of use
typedef Curve<2, float> Curve2f;
typedef Curve<2, double> Curve2d;
typedef Curve<3, float> Curve3f;
typedef Curve<3, double> Curve3d;
typedef RationalCurve<2, float> RationalCurve2f;
typedef RationalCurve<2, double> RationalCurve2d;
typedef RationalCurve<3, float> RationalCurve3f;
typedef RationalCurve<3, double> RationalCurve3d;

} // namespace tinynurbs

#endif // TINYNURBS_CURVE_H
