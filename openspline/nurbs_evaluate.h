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
bool isValidRelation(unsigned int degree, size_t nKnots, size_t nCtrlPts);

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

	// Find span and corresponding non-zero basis functions
	int span = findSpan(degree, knots, u);
	std::vector<double> N;
	bsplineBasis(degree, span, knots, u, N);

	// Compute point
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

/**
Evaluate derivatives of a non-rational NURBS curve
@param[in] u Parameter to evaluate the derivatives at.
@param[in] knots Knot vector of the curve.
@param[in] controlPoints Control points of the curve.
@param[in] nDers Number of times to derivate.
@param[in, out] curveDers Derivatives of the curve.
E.g. curveDers[n] is the nth derivative at u, where 0 <= n <= nDers.
*/
template <int dim, typename T>
void nurbsCurveDerivatives(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec<dim, T>> &controlPoints,
	int nDers, std::vector<glm::vec<dim, T>> &curveDers) {

	typedef glm::vec<dim, T> tvecn;
	using std::vector;

	curveDers.clear();
	curveDers.resize(nDers + 1);

	// Assign higher order derivatives to zero
	for (int k = degree + 1; k <= nDers; k++) {
		curveDers[k] = tvecn(0.0);
	}
	
	// Find the span and corresponding non-zero basis functions & derivatives
	int span = findSpan(degree, knots, u);
	vector<vector<double>> ders;
	bsplineDerBasis(degree, span, knots, u, nDers, ders);

	// Compute first nDers derivatives
	int du = nDers < degree ? nDers : degree;
	for (int k = 0; k <= du; k++) {
		curveDers[k] = tvecn(0.0);
		for (int j = 0; j <= degree; j++) {
			curveDers[k] += static_cast<T>(ders[k][j]) *
				controlPoints[span - degree + j];
		}
	}
}

/**
Evaluate derivatives of a rational NURBS curve
@param[in] u Parameter to evaluate the derivatives at.
@param[in] knots Knot vector of the curve.
@param[in] controlPoints Control points of the curve.
@param[in] weights Weights corresponding to each control point.
@param[in] nDers Number of times to derivate.
@param[in, out] curveDers Derivatives of the curve.
E.g. curveDers[n] is the nth derivative at u, where n is between 0 and nDers-1.
*/
template <int dim, typename T>
void nurbsRationalCurveDerivatives(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec<dim, T>> &controlPoints,
	const std::vector<T> weights, int nDers,
	std::vector<glm::vec<dim, T>> &curveDers) {

	typedef glm::vec<dim, T> tvecn;
	typedef glm::vec<dim + 1, T> tvecnp1;
	using std::vector;

	// Compute homogenous coordinates of control points
	vector<tvecnp1> Cw;
	Cw.reserve(controlPoints.size());
	for (int i = 0; i < controlPoints.size(); i++) {
		Cw.push_back(tvecnp1(
			util::cartesianToHomogenous(controlPoints[i], weights[i])
		));
	}

	// Derivatives of Cw
	vector<tvecnp1> Cwders;
	nurbsCurveDerivatives(u, degree, knots, Cw, nDers, Cwders);

	// Split Cwders into coordinates and weights
	vector<tvecn> Aders;
	vector<T> wders;
	for (const auto &val : Cwders) {
		tvecn Aderi;
		for (int i = 0; i < dim - 1; i++) {
			Aderi[i] = val[i];
		}
		Aders.push_back(Aderi);
		wders.push_back(val[dim - 1]);
	}

	// Compute rational derivatives
	for (int k = 0; k <= nDers; k++) {
		tvecn v = Aders[k];
		for (int i = 1; i <= k; i++) {
			v -= static_cast<T>(util::binomial(k, i)) * wders[i] * curveDers[k - i];
		}
		curveDers[k] = v / wders[0];
	}
}

/**
Evaluate point on a nonrational NURBS surface
@param[in] u Parameter to evaluate the surface at.
@param[in] v Parameter to evaluate the surface at.
@param[in] degreeU Degree of the given surface in u-direction.
@param[in] degreeV Degree of the given surface.
@param[in] knots Knot vector of the curve.
@param[in] controlPoints Control points of the curve.
@param[in, out] point Resulting point on the curve at parameter u.
*/
template <int dim, typename T>
void nurbsSurfacePoint(double u, double v, uint8_t degreeU, uint8_t degreeV,
	const std::vector<double> &knotsU, const std::vector<double> &knotsV,
	const std::vector<std::vector<glm::vec<dim, T>>> &controlPoints,
	glm::vec<dim, T> &point) {

	// Initialize result to 0s
	for (int i = 0; i < dim; i++) {
		point[i] = static_cast<T>(0.0);
	}

	// Find span and non-zero basis functions
	int spanU = findSpan(degreeU, knotsU, u);
	int spanV = findSpan(degreeV, knotsV, v);
	std::vector<double> Nu, Nv;
	bsplineBasis(degreeU, spanU, knotsU, u, Nu);
	bsplineBasis(degreeV, spanV, knotsV, v, Nv);

	for (int l = 0; l <= degreeV; l++) {
		glm::vec<dim, T> temp(0.0);
		for (int k = 0; k <= degreeU; k++) {
			temp += static_cast<T>(Nu[k]) *
				controlPoints[spanU - degreeU + k][spanV - degreeV + l];
		}

		point += static_cast<T>(Nv[l]) * temp;
	}
}


/**
Evaluate point on a non-rational NURBS surface
@param[in] u Parameter to evaluate the surface at.
@param[in] v Parameter to evaluate the surface at.
@param[in] degreeU Degree of the given surface in u-direction.
@param[in] degreeV Degree of the given surface.
@param[in] knots Knot vector of the curve.
@param[in] controlPoints Control points of the curve.
@param[in, out] point Resulting point on the curve at parameter u.
*/
template <int dim, typename T>
void nurbsRationalSurfacePoint(double u, double v, uint8_t degreeU, uint8_t degreeV,
	const std::vector<double> &knotsU, const std::vector<double> &knotsV,
	const std::vector<std::vector<glm::vec<dim, T>>> &controlPoints,
	const std::vector<std::vector<T>> &weights,
	glm::vec<dim, T> &point) {

}

/**
Evaluate derivatives on a non-rational NURBS surface
@param[in] u Parameter to evaluate the surface at.
@param[in] v Parameter to evaluate the surface at.
@param[in] degreeU Degree of the given surface in u-direction.
@param[in] degreeV Degree of the given surface.
@param[in] knots Knot vector of the curve.
@param[in] controlPoints Control points of the curve.
@param[in, out] surfDers Resulting point on the curve at parameter u.
*/
template <int dim, typename T>
void nurbsSurfaceDerivatives(double u, double v,
	uint8_t degreeU, uint8_t degreeV,
	const std::vector<double> &knotsU, const std::vector<double> &knotsV,
	const std::vector<std::vector<glm::vec<dim, T>>> &controlPoints,
	int nDers, std::vector<std::vector<glm::vec<dim, T>>> &surfDers) {

	surfDers.clear();
	surfDers.resize(nDers + 1);
	for (auto &vec : surfDers) {
		vec.resize(nDers + 1);
	}

	// Set higher order derivatives to 0
	for (int k = degreeU + 1; k <= nDers; k++) {
		for (int l = degreeV + 1; l <= nDers; l++) {
			surfDers[k][l] = glm::vec<dim, T>(0.0);
		}
	}

	// Find span and basis function derivatives
	int spanU = findSpan(degreeU, knotsU, u);
	int spanV = findSpan(degreeV, knotsV, v);
	std::vector<std::vector<double>> dersU, dersV;
	bsplineDerBasis(degreeU, spanU, knotsU, u, nDers, dersU);
	bsplineDerBasis(degreeV, spanV, knotsV, v, nDers, dersV);

	// Number of non-zero derivatives is <= degree
	int du = nDers < degreeU ? nDers : degreeU;
	int dv = nDers < degreeV ? nDers : degreeV;
	std::vector<glm::vec<dim, T>> temp;
	temp.resize(degreeV + 1);

	// Compute derivatives
	for (int k = 0; k <= du; k++) {
		for (int s = 0; s <= degreeV; s++) {
			temp[s] = glm::vec<dim, T>(0.0);
			for (int r = 0; r <= degreeU; r++) {
				temp[s] += static_cast<T>(dersU[k][r]) * 
					controlPoints[spanU - degreeU + r][spanV - degreeV + s];
			}
		}

		int nk = nDers - k;
		int dd = nk < dv ? nk : dv;

		for (int l = 0; l <= dd; l++) {
			surfDers[k][l] = glm::vec<dim, T>(0.0);

			for (int s = 0; s <= degreeV; s++) {
				surfDers[k][l] += static_cast<T>(dersV[l][s]) * temp[s];
			}
		}
	}
}


template <int dim, typename T>
void nurbsRationalSurfaceDerivatives(double u, double v,
	uint8_t degreeU, uint8_t degreeV,
	const std::vector<double> &knotsU, const std::vector<double> &knotsV,
	const std::vector<std::vector<glm::vec<dim, T>>> &controlPoints,
	const std::vector<std::vector<T>> &weights,
	int nDers, std::vector<std::vector<glm::vec<dim, T>>> &surfDers) {

}

} // namespace ospl
