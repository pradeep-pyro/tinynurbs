#include "nurbs_evaluate.h"
#include "array2.h"

// Explicit template function instantiations for glm vec2, vec3, dvec2
// and dvec3 types; these are the only relevant types.
template void ospl::nurbs::curvePoint<2, float>(double u, uint8_t degree, 
	const std::vector<double> &knots,
	const std::vector<glm::vec2> &controlPoints, glm::vec2 &point);
template void ospl::nurbs::curvePoint<2, double>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec2> &controlPoints, glm::dvec2 &point);
template void ospl::nurbs::curvePoint<3, float>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec3> &controlPoints, glm::vec3 &point);
template void ospl::nurbs::curvePoint<3, double>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec3> &controlPoints, glm::dvec3 &point);

template void ospl::nurbs::rationalCurvePoint<2, float>(double u, uint8_t degree,
	const std::vector<double> &knots, 
	const std::vector<glm::vec2> &controlPoints,
	const std::vector<float> &weights, glm::vec2 &point);
template void ospl::nurbs::rationalCurvePoint<2, double>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec2> &controlPoints,
	const std::vector<double> &weights, glm::dvec2 &point);
template void ospl::nurbs::rationalCurvePoint<3, float>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec3> &controlPoints,
	const std::vector<float> &weights, glm::vec3 &point);
template void ospl::nurbs::rationalCurvePoint<3, double>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec3> &controlPoints,
	const std::vector<double> &weights, glm::dvec3 &point);

template void ospl::nurbs::curveDerivatives<2, float>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec2> &controlPoints,
	int nDers, std::vector<glm::vec2> &ders);
template void ospl::nurbs::curveDerivatives<2, double>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec2> &controlPoints,
	int nDers, std::vector<glm::dvec2> &ders);
template void ospl::nurbs::curveDerivatives<3, float>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::vec3> &controlPoints,
	int nDers, std::vector<glm::vec3> &ders);
template void ospl::nurbs::curveDerivatives<3, double>(double u, uint8_t degree,
	const std::vector<double> &knots,
	const std::vector<glm::dvec3> &controlPoints,
	int nDers, std::vector<glm::dvec3> &ders);