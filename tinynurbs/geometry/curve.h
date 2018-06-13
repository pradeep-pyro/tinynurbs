/*
@file nurbstk/curve.h
@author Pradeep Kumar Jayaraman <pradeep.pyro@gmail.com>

The Curve class represents a non-uniform polynomial B-spline curve, while the RationalCurve class
represents a non-uniform rational B-spline (NURBS) curve.

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE.txt file.
*/

#pragma once
#include <vector>
#include <exception>
#include <stdexcept>
#include "glm/glm.hpp"

namespace tinynurbs {

/**
Checks if the relation between degree, number of knots, and
number of control points is valid
@param degree Degree of the NURBS curve
@param num_knots Number of knot values
@param num_ctrl_pts Number of control points
@return Whether the relationship is valid
*/
bool isValidRelation(unsigned int degree, size_t num_knots, size_t num_ctrl_pts);

/**
@brief Struct for holding a polynomial B-spline curve
@tparam dim Dimension of the curve (2 or 3)
@tparam T Data type of control points and knots (float or double)
*/
template <int dim, typename T>
struct Curve {
    unsigned int degree;
    std::vector<T> knots;
    std::vector<glm::vec<dim, T>> control_points;

    Curve(unsigned int degree, const std::vector<T> &knots,
          const std::vector<glm::vec<dim, T>> &control_points)
        : degree(degree), knots(knots), control_points(control_points) {
    }
};

/**
@brief Struct for holding a rational B-spline curve
@tparam dim Dimension of the curve (2 or 3)
@tparam T Data type of control points and knots (float or double)
*/
template <int dim, typename T>
struct RationalCurve : public Curve<dim, T> {
    std::vector<T> weights;

    RationalCurve(unsigned int degree, const std::vector<T> &knots,
                  const std::vector<glm::vec<dim, T>> &control_points,
                  const std::vector<T> weights)
        : weights(weights), Curve<dim, T>(degree, knots, control_points) {
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
