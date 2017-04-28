#pragma once

#include <vector>

namespace ospl {
namespace nurbs {

int findSpan(int degree, const std::vector<double> &knots, double u);

double basisFunction(int i, int deg, const std::vector<double> &U, double u);

void basisFunctions(int deg, int span, const std::vector<double> &knots,
	double u, std::vector<double> &N);

void derivativeBasisFunctions(int deg, int span, const std::vector<double> &knots,
	double u, int nDers, std::vector<std::vector<double>> &Nk);

} // namespace nurbs
} // namespace ospl