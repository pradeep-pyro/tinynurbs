/*
@file openspline/nurbs_evaluate.h
@author Pradeep Kumar Jayaraman <pradeep.pyro@gmail.com>

Core functionality for evaluating points and derivatives on NURBS curves and
surfaces

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE.txt file.
*/

#pragma once

#include <vector>

#include "glm/glm.hpp"

namespace {

bool close(double a, double b, double eps = std::numeric_limits<double>::epsilon()) {
	return (std::abs(a - b) < eps) ? true : false;
}

} // namespace


namespace ospl {
namespace nurbs {

int findSpan(int degree, const std::vector<double> &knots, double u);

double basisFunction(int i, int p, const std::vector<double> &U, double u);

void basisFunctions(int deg, int span, const std::vector<double> &knots,
	double u, std::vector<double> &N);

/**
Evaluate point on a nonrational NURBS curve
@param[in] u Parameter to evaluate the curve at.
@param[in] degree Degree of the given curve.
@param[in] knots Knot vector of the curve.
@param[in] controlPoints Control points of the curve.
@param[in, out] point Resulting point on the curve at parameter u.
*/
template <int dim, typename T>
void curvePoint(double u, int degree, const std::vector<double> &knots,
	const std::vector<glm::vec<dim, T>> &controlPoints, glm::vec<dim, T> &point) {
	// Initialize result to 0s
	for (int i = 0; i < dim; i++) {
		point[i] = static_cast<T>(0.0);
	}

	// Find span and non-zero basis functions
	int span = findSpan(degree, knots, u);
	std::vector<double> N;
	basisFunctions(degree, span, knots, u, N);

	// Compute point: \sum_{j = 0}^{degree} N_j * C_j
	for (int j = 0; j <= degree; j++) {
		point += static_cast<T>(N[j]) * controlPoints[span - degree + j];
	}
}

/**
Convert an nd point in homogenous coordinates to an (n-1)d point in cartesian
coordinates by perspective division
@param[in] pt Point in homogenous coordinates
@return Input point in cartesian coordinates
*/
template<int nd, typename T>
glm::vec<nd - 1, T> homogenousToCartesian(glm::vec<nd, T> pt) {
	return glm::vec<nd - 1, T>(pt / pt[pt.length()-1]);
}

/**
Convert an nd point in cartesian coordinates to an (n+1)d point in homogenous
coordinates
@param[in] pt Point in cartesian coordinates
@param[in] w Weight
@return Input point in homogenous coordinates
*/
template<int nd, typename T>
glm::vec<nd + 1, T> cartesianToHomogenous(glm::vec<nd, T> pt, T w) {
	return glm::vec<nd + 1, T>(pt * w, w);
}

/**
Evaluate point on a rational NURBS curve
@param[in] u Parameter to evaluate the curve at.
@param[in] knots Knot vector of the curve.
@param[in] controlPoints Control points of the curve.
@param[in] weights Weights corresponding to each control point.
@param[in, out] point Resulting point on the curve.
*/
template <int dim, typename T>
void rationalCurvePoint(double u, int degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec<dim, T>> &controlPoints,
	const std::vector<T> &weights, glm::vec<dim, T> &point) {

	typedef glm::vec<dim + 1, T> tvecnp1;

	std::vector<tvecnp1> Cw;
	Cw.reserve(controlPoints.size());
	for (int i = 0; i < controlPoints.size(); i++) {
		Cw.push_back(tvecnp1(
					cartesianToHomogenous(controlPoints[i], weights[i])
			        ));
	}

	tvecnp1 pointw;
	curvePoint(u, degree, knots, Cw, pointw);
	point = homogenousToCartesian(pointw);
}

} // namespace nurbs
} // namespace ospl
