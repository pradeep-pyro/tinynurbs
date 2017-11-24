#include "evaluate.h"
#include "util.h"

bool nurbstk::isValidRelation(unsigned int degree, size_t nKnots, size_t nCtrlPts) {
	return nKnots == nCtrlPts + degree + 1;
}

// Explicit template function instantiations for glm vec2, vec3, dvec2
// and dvec3 types; these are the only relevant types.
template void nurbstk::nurbsCurvePoint<2, float>(double u, uint8_t degree, 
	const std::vector<double> &knots,
	const std::vector<glm::vec2> &controlPoints, glm::vec2 &point);
template void nurbstk::nurbsCurvePoint<2, double>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec2> &controlPoints, glm::dvec2 &point);
template void nurbstk::nurbsCurvePoint<3, float>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec3> &controlPoints, glm::vec3 &point);
template void nurbstk::nurbsCurvePoint<3, double>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec3> &controlPoints, glm::dvec3 &point);

template void nurbstk::nurbsRationalCurvePoint<2, float>(double u, uint8_t degree,
	const std::vector<double> &knots, 
	const std::vector<glm::vec2> &controlPoints,
	const std::vector<float> &weights, glm::vec2 &point);
template void nurbstk::nurbsRationalCurvePoint<2, double>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec2> &controlPoints,
	const std::vector<double> &weights, glm::dvec2 &point);
template void nurbstk::nurbsRationalCurvePoint<3, float>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec3> &controlPoints,
	const std::vector<float> &weights, glm::vec3 &point);
template void nurbstk::nurbsRationalCurvePoint<3, double>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec3> &controlPoints,
	const std::vector<double> &weights, glm::dvec3 &point);

template void nurbstk::nurbsCurveDerivatives<2, float>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec2> &controlPoints,
	int nDers, std::vector<glm::vec2> &ders);
template void nurbstk::nurbsCurveDerivatives<2, double>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec2> &controlPoints,
	int nDers, std::vector<glm::dvec2> &ders);
template void nurbstk::nurbsCurveDerivatives<3, float>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec3> &controlPoints,
	int nDers, std::vector<glm::vec3> &ders);
template void nurbstk::nurbsCurveDerivatives<3, double>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec3> &controlPoints,
	int nDers, std::vector<glm::dvec3> &ders);

template void nurbstk::nurbsRationalCurveDerivatives<2, float>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec2> &controlPoints,
	const std::vector<float> weights,
	int nDers, std::vector<glm::vec2> &curveDers);
template void nurbstk::nurbsRationalCurveDerivatives<2, double>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec2> &controlPoints,
	const std::vector<double> weights,
	int nDers, std::vector<glm::dvec2> &curveDers);
template void nurbstk::nurbsRationalCurveDerivatives<3, float>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec3> &controlPoints,
	const std::vector<float> weights,
	int nDers, std::vector<glm::vec3> &curveDers);
template void nurbstk::nurbsRationalCurveDerivatives<3, double>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec3> &controlPoints,
	const std::vector<double> weights,
	int nDers, std::vector<glm::dvec3> &curveDers);

template void nurbstk::nurbsSurfacePoint<3, float>(double u, double v,
	uint8_t degreeU, uint8_t degreeV,
	const std::vector<double> &knotsU, const std::vector<double> &knotsV,
	const std::vector<std::vector<glm::vec3>> &controlPoints,
	glm::vec3 &point);
template void nurbstk::nurbsSurfacePoint<3, double>(double u, double v,
	uint8_t degreeU, uint8_t degreeV,
	const std::vector<double> &knotsU, const std::vector<double> &knotsV,
	const std::vector<std::vector<glm::dvec3>> &controlPoints,
	glm::dvec3 &point);