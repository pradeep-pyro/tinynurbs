#pragma once

#include <vector>

namespace ospl {

int findSpan(int degree, const std::vector<double> &knots, double u);

double bsplineOneBasis(int i, int deg, const std::vector<double> &U, double u);

void bsplineBasis(int deg, int span, const std::vector<double> &knots,
	double u, std::vector<double> &N);

void bsplineDerBasis(int deg, int span, const std::vector<double> &knots,
	double u, int nDers, std::vector<std::vector<double>> &Nk);

} // namespace ospl