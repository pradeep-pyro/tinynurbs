/*
@file tinynurbs/core/check.h
@author Pradeep Kumar Jayaraman <pradeep.pyro@gmail.com>

Functions to refine curves and surfaces

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE.txt file.
*/

#ifndef REFINE_H
#define REFINE_H

#include <vector>
#include "glm/glm.hpp"
#include "../core/check.h"
#include "../util/util.h"
#include "../geometry/curve.h"

namespace tinynurbs {

template <int dim, typename T>
void curveKnotRefine(unsigned int deg, std::vector<T> &knots, std::vector<glm::vec<dim, T>> &cp,
                     T u, unsigned int r=1) {
    int k = findSpan(deg, knots, u);
    unsigned int s = knotMultiplicity(knots, k);
    // Insert new knots between span and (span + 1)
    std::vector<T> new_knots;
    new_knots.resize(knots.size() + r);
    for (int i = 0; i < k + 1; ++i) {
        new_knots[i] = knots[i];
    }
    for (int i = 1; i < r + 1; ++i) {
        new_knots[k + i] = u;
    }
    for (int i = k + 1; i < knots.size(); ++i) {
        new_knots[i + r] = knots[i];
    }
    // Copy unaffected control points
    std::vector<glm::vec<dim, T>> new_cp;
    new_cp.resize(cp.size() + r);
    for (int i = 0; i < k - deg + 1; ++i) {
        new_cp[i] = cp[i];
    }
    for (int i = k - s; i < cp.size(); ++i) {
        new_cp[i + r] = cp[i];
    }
    // Copy affected control points
    std::vector<glm::vec<dim, T>> tmp;
    tmp.resize(deg - s + 1);
    for (int i = 0; i < deg - s + 1; ++i) {
        tmp[i] = cp[k - deg + i];
    }
    // Modify affected control points
    for (int j = 1; j < r + 1; ++j) {
        int L = k - deg + j;
        for (int i = 0; i < deg - j - s + 1; ++i) {
            T a = (u - knots[L + i]) / (knots[i + k + 1] - knots[L + i]);
            tmp[i] = (1 - a) * tmp[i] + a * tmp[i + 1];
        }
        new_cp[L] = tmp[0];
        new_cp[k + r - j - s] = tmp[deg - j - s];
    }
    int L = k - deg + r;
    for (int i = L + 1; i < k - s; ++i) {
        new_cp[i] = tmp[i - L];
    }
    // Update original control points and knots
    cp = new_cp;
    knots = new_knots;
}

template <int dim, typename T>
void curveKnotRefine(Curve<dim, T> &crv, T u, unsigned int r=1) {
    curveKnotRefine(crv.degree, crv.knots, crv.control_points, u, r);
}

template <int dim, typename T>
void curveKnotRefine(RationalCurve<dim, T> &crv, T u, unsigned int r=1) {
    std::vector<glm::vec<dim + 1, T>> Cw;
    Cw.reserve(crv.control_points.size());
    for (int i = 0; i < crv.control_points.size(); ++i) {
        Cw.push_back(util::cartesianToHomogenous(crv.control_points[i], crv.weights[i]));
    }
    curveKnotRefine(crv.degree, crv.knots, Cw, u, r);
    std::vector<glm::vec<dim, T>> new_cp, new_w;
    new_cp.reserve(Cw.size());
    new_w.reserve(Cw.size());
    for (int i = 0; i < crv.control_points.size(); ++i) {
        new_cp.push_back(util::homogenousToCartesian(Cw[i]));
        new_w.push_back(Cw[dim - 1]);
    }
    crv.control_points = new_cp;
    crv.weights = new_w;
}

} // namespace tinynurbs

#endif // REFINE_H
