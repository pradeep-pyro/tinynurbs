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

#include "bspline_basis.h"
#include "util.h"

namespace ospl {

/**
Checks if the relation between degree, number of knots, and
number of control points is valid
@param degree Degree of the NURBS curve
@param nKnots Number of knot values
@param nCtrlPts Number of control points
@return Whether the relationship is valid
*/
bool isValidRelation(int degree, int nKnots, int nCtrlPts) {
	if (nKnots != nCtrlPts + degree + 1) {
		return false;
	}
	return true;
}


/**
Evaluate point on a nonrational NURBS curve
@param[in] u Parameter to evaluate the curve at.
@param[in] degree Degree of the given curve.
@param[in] knots Knot vector of the curve.
@param[in] controlPoints Control points of the curve.
@param[in, out] point Resulting point on the curve at parameter u.
*/
template <int dim, typename T>
void nurbsCurvePoint(double u, uint8_t degree, const std::vector<double> &knots,
	const std::vector<glm::vec<dim, T>> &controlPoints, glm::vec<dim, T> &point) {
	// Initialize result to 0s
	for (int i = 0; i < dim; i++) {
		point[i] = static_cast<T>(0.0);
	}

	// Find span and non-zero basis functions
	int span = findSpan(degree, knots, u);
	std::vector<double> N;
	bsplineBasis(degree, span, knots, u, N);

	// Compute point: \sum_{j = 0}^{degree} N_j * C_j
	for (int j = 0; j <= degree; j++) {
		point += static_cast<T>(N[j]) * controlPoints[span - degree + j];
	}
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
void nurbsRationalCurvePoint(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec<dim, T>> &controlPoints,
	const std::vector<T> &weights, glm::vec<dim, T> &point) {

	typedef glm::vec<dim + 1, T> tvecnp1;

	// Compute homogenous coordinates of control points
	std::vector<tvecnp1> Cw;
	Cw.reserve(controlPoints.size());
	for (int i = 0; i < controlPoints.size(); i++) {
		Cw.push_back(tvecnp1(
			util::cartesianToHomogenous(controlPoints[i], weights[i])
		));
	}

	// Compute point using homogenous coordinates
	tvecnp1 pointw;
	nurbsCurvePoint(u, degree, knots, Cw, pointw);

	// Convert back to cartesian coordinates
	point = util::homogenousToCartesian(pointw);
}

template <int dim, typename T>
void nurbsCurveDerivatives(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec<dim, T>> &controlPoints,
	int nDers, std::vector<glm::vec<dim, T>> &curveDers) {

	curveDers.clear();
	curveDers.resize(nDers + 1);
	for (int k = degree + 1; k <= nDers; k++) {
		curveDers[k] = glm::vec<dim, T>(0.0);
	}
	
	int span = findSpan(degree, knots, u);
	std::vector<std::vector<double>> nders;
	bsplineDerBasis(degree, span, knots, u, nDers, nders);
	int du = nDers < degree ? nDers : degree;
	for (int k = 0; k <= du; k++) {
		curveDers[k] = glm::vec<dim, T>(0.0);
		for (int j = 0; j <= degree; j++) {
			curveDers[k] += static_cast<T>(nders[k][j]) *
				controlPoints[span - degree + j];
		}
	}
}

} // namespace ospl
