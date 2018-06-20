/*
@file nurbstk/nurbs_knots.h
@author Pradeep Kumar Jayaraman <pradeep.pyro@gmail.com>

Helper functions for creating and modifying knot vectors.

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE.txt file.
*/

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
