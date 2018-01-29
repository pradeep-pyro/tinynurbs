#include "knots.h"
/*
#include <iostream>
#include <algorithm>

using namespace std;

bool makeUniformKnotVector(unsigned int degree, size_t nCtrlPts,
    std::vector<double> &knots) {
    // Compute number of knots
    if (nCtrlPts < degree + 1) {
        return false;
    }
    int nKnots = nCtrlPts + degree + 1;

    // Interior knots
    double step = 1.0 / (nKnots - 1);
    knots.clear();
    knots.reserve(nKnots);
    for (double u = 0.0; u <= 1.0; u += step) {
        knots.push_back(u);
    }

    return true;
}

bool makeClampedUniformKnotVector(unsigned int degree, size_t nCtrlPts,
    std::vector<double> &knots) {
    // Compute number of knots
    if (nCtrlPts < degree + 1) {
        return false;
    }
    int nKnots = nCtrlPts + degree + 1;
    int nIntKnots = nKnots - 2 * degree;

    knots.clear();
    knots.reserve(nKnots);

    // Clamp left side
    for (int i = 0; i < degree; i++) {
        knots.push_back(0.0);
    }

    // Interior knots
    if (nIntKnots > 0) {
        double step = 1.0 / (nIntKnots - 1);
        double u = 0.0;
        for (int i = 0; i < nIntKnots; i++) {
            knots.push_back(static_cast<double>(i) * step);
        }
    }

    // Clamp right side
    for (int i = 0; i < degree; i++) {
        knots.push_back(1.0);
    }

    return true;
}

void clampKnotVector(unsigned int degree, std::vector<double>& knots) {
    clampKnotVectorLeft(degree, knots);
    clampKnotVectorRight(degree, knots);
}

void clampKnotVectorLeft(unsigned int degree, std::vector<double>& knots) {
    double start = knots[degree];
    for (int i = 0; i < degree; i++) {
        knots[i] = start;
    }
}

void clampKnotVectorRight(unsigned int degree, std::vector<double> &knots) {
    double end = knots[knots.size() - degree - 1];
    for (int i = 0; i < degree; i++) {
        knots[knots.size() - 1 - i] = end;
    }
}

bool isKnotVectorMonotonic(const std::vector<double> &knots) {
    return std::is_sorted(knots.begin(), knots.end());
}

bool isKnotVectorClosed(unsigned int degree, const std::vector<double> &knots) {
    for (int i = 0; i < degree + 2; i++) {
        if (abs(knots[i] - knots[knots.size() - degree - 2 + i])
            > numeric_limits<double>::epsilon()) {
            return false;
        }
    }
    return true;
}
*/
