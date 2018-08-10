/*
@file tinynurbs/modify/nurbs_knots.h
@author Pradeep Kumar Jayaraman <pradeep.pyro@gmail.com>

Functions to modify knot vectors

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
bool makeUniformKnotVector(unsigned int degree, size_t num_ctrl_pts,
                           std::vector<T> &knots) {
    // Compute number of knots
    if (num_ctrl_pts < degree + 1) {
        return false;
    }
    int num_knots = num_ctrl_pts + degree + 1;

    // Interior knots
    T step = 1.0 / (num_knots - 1);
    knots.clear();
    knots.reserve(num_knots);
    for (T u = 0.0; u <= 1.0; u += step) {
        knots.push_back(u);
    }

    return true;
}

template <typename T>
bool makeClampedUniformKnotVector(unsigned int degree, size_t num_ctrl_pts,
                                  std::vector<T> &knots) {
    // Compute number of knots
    if (num_ctrl_pts < degree + 1) {
        return false;
    }
    int num_knots = num_ctrl_pts + degree + 1;
    int num_int_knots = num_knots - 2 * degree;

    knots.clear();
    knots.reserve(num_knots);

    // Clamp left side
    for (int i = 0; i < degree; i++) {
        knots.push_back(0.0);
    }

    // Interior knots
    if (num_int_knots > 0) {
        double step = 1.0 / (num_int_knots - 1);
        for (int i = 0; i < num_int_knots; i++) {
            knots.push_back(static_cast<double>(i) * step);
        }
    }

    // Clamp right side
    for (int i = 0; i < degree; i++) {
        knots.push_back(1.0);
    }

    return true;
}

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
