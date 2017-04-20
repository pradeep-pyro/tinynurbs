#include "nurbs_evaluate.h"


namespace ospl {
namespace nurbs {

/**
Find the span of the given parameter in the knot vector.
@param[in] degree Degree of the curve.
@param[in] knots Knot vector of the curve.
@param[in] u Parameter value.
@return Span index into the knot vector such that (span - 1) < u <= span
*/
int findSpan(int degree, const std::vector<double> &knots, double u) {
	int n = knots.size() - degree - 2; // index of last control point

	// For u that is equal to last knot value
	if (close(u, knots[n + 1])) {
		return n;
	}

	// For values of u that lies outside the domain
	if (u > knots[n + 1]) {
		return n;
	}
	if (u < knots[degree]) {
		return degree;
	}

	// Binary search
	// TODO: Replace this with std::lower_bound
	int low = degree;
	int high = n + 1;
	int mid = (int)std::floor((low + high) / 2.0);
	while (u < knots[mid] || u >= knots[mid + 1]) {
		if (u < knots[mid]) {
			high = mid;
		}
		else {
			low = mid;
		}
		mid = (low + high) / 2;
	}
	return mid;
}

/**
Compute a single B-spline basis function
*/
// 
double basisFunction(int i, int p, const std::vector<double> &U, double u) {
	int m = U.size() - 1;
	// Special case
	if ((i == 0 && close(u, U[0])) || (i == m - p - 1 && close(u, U[m]))) {
		return 1.0;
	}
	// Local Property
	if (u < U[i] || u >= U[i + p + 1]) {
		return 0.0;
	}
	// Initialize zeroth-degree functions
	std::vector<double> N;
	N.resize(p + 1);
	for (int j = 0; j <= p; j++) {
		N[j] = (u >= U[i + j] && u < U[i + j + 1]) ? 1.0 : 0.0;
	}
	// Compute triangular table
	for (int k = 1; k <= p; k++) {
		double saved = (close(N[0], 0.0)) ? 0.0
			: ((u - U[i])*N[0]) / (U[i + k] - U[i]);
		for (int j = 0; j < p - k + 1; j++) {
			double Uleft = U[i + j + 1];
			double Uright = U[i + j + k + 1];
			if (close(N[j + 1], 0.0)) {
				N[j] = saved;
				saved = 0.0;
			}
			else {
				double temp = N[j + 1] / (Uright - Uleft);
				N[j] = saved + (Uright - u)*temp;
				saved = (u - Uleft)*temp;
			}
		}
	}
	return N[0];
}

/**
// Compute all non-zero B-spline basis functions
@param[in] deg Degree of the basis function.
@param[in] span Index obtained from findSpan() corresponding the u and knots.
@param[in] knots Knot vector corresponding to the basis functions.
@param[in] u Parameter to evaluate the basis functions at.
@param[in, out] N Values of (deg+1) non-zero basis functions.
*/
void basisFunctions(int deg, int span, const std::vector<double> &knots, double u,
	std::vector<double> &N) {
	N.clear();
	N.resize(deg + 1, 0.0);
	std::vector<double> left, right;
	left.resize(deg + 1, 0.0);
	right.resize(deg + 1, 0.0);
	double saved = 0.0, temp = 0.0;

	N[0] = 1.0;

	for (int j = 1; j <= deg; j++) {
		left[j] = (u - knots[span + 1 - j]);
		right[j] = knots[span + j] - u;
		saved = 0.0;
		for (int r = 0; r < j; r++) {
			temp = N[r] / (right[r + 1] + left[j - r]);
			N[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		N[j] = saved;
	}
}

} // namespace nurbs
} // namespace ospl

// Explicit template function instantiations for glm vec2, vec3, dvec2
// and dvec3 types; these are the only relevant types.
template void ospl::nurbs::curvePoint<2, float>(double u, int degree, 
	const std::vector<double> &knots,
	const std::vector<glm::vec2> &controlPoints, glm::vec2 &point);
template void ospl::nurbs::curvePoint<2, double>(double u, int degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec2> &controlPoints, glm::dvec2 &point);
template void ospl::nurbs::curvePoint<3, float>(double u, int degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec3> &controlPoints, glm::vec3 &point);
template void ospl::nurbs::curvePoint<3, double>(double u, int degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec3> &controlPoints, glm::dvec3 &point);

template void ospl::nurbs::rationalCurvePoint<2, float>(double u, int degree,
	const std::vector<double> &knots, 
	const std::vector<glm::vec2> &controlPoints,
	const std::vector<float> &w, glm::vec2 &point);
template void ospl::nurbs::rationalCurvePoint<2, double>(double u, int degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec2> &controlPoints,
	const std::vector<double> &w, glm::dvec2 &point);
template void ospl::nurbs::rationalCurvePoint<3, float>(double u, int degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec3> &controlPoints,
	const std::vector<float> &w, glm::vec3 &point);
template void ospl::nurbs::rationalCurvePoint<3, double>(double u, int degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec3> &controlPoints,
	const std::vector<double> &w, glm::dvec3 &point);