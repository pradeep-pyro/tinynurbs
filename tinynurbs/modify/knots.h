#ifndef KNOTS_H
#define KNOTS_H

#include <algorithm>
#include <vector>
#include "glm/glm.hpp"

namespace tinynurbs {

template <typename T>
void clampKnotVector(unsigned int degree, std::vector<T> &knots) {
    clampKnotVectorLeft(degree, knots);
    clampKnotVectorRight(degree, knots);
}

template <typename T>
void clampKnotVectorLeft(unsigned int degree, std::vector<T> &knots) {
    T start = knots[degree];
    for (int i = 0; i < degree; i++) {
        knots[i] = start;
    }
}

template <typename T>
void clampKnotVectorRight(unsigned int degree, std::vector<T> &knots) {
    T end = knots[knots.size() - degree - 1];
    for (int i = 0; i < degree; i++) {
        knots[knots.size() - 1 - i] = end;
    }
}

} // namespace tinynurbs

#endif
