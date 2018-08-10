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

template <typename T, int n>
void curveInsertKnot(T u, unsigned int deg, std::vector<glm::vec<T, n>> &knots, std::vector<glm::vec<T, n>> &cp) {
    int span = findSpan(u, deg, knots);
    // Insert new control point and knot between span and (span + 1)
    knots.insert(span, u);
    cp.insert(span, glm::vec<T, n>(0.0));
    // Compute new positions for affected control points
    for (int i = span - deg + 1; i <= span; ++i) {
        T a = (u - knots[i]) / (knots[i + deg] - knots[i]);
        cp[i] = (1 - a) * cp[i - 1] + a * cp[i];
    }
}

#endif // REFINE_H
