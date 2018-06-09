#ifndef REFINE_H
#define REFINE_H

#include <vector>
#include "glm/glm.hpp"

template <typename T, int n>
void curveInsertKnot(T u, unsigned int deg, std::vector<glm::vec<T, n>> &knots, std::vector<glm::vec<T, n>> &cp) {
    int span = findSpan(u, deg, knots);
    // Insert new control point and knot between span and (span + 1)
    knots.insert(span, u);
    cp.insert(span, tvecn(0.0));
    // Compute new positions for affected control points
    for (int i = span - deg + 1; i <= span; ++i) {
        T a = (u - knots[i]) / (knots[i + deg] - knots[i]);
        cp[i] = (1 - a) * cp[i - 1] + a * cp[i];
    }
}

template <typename T, int n>
void surfaceInsertKnotU(T u, unsigned int deg, std::vector<glm::vec<T, n>> &knots, std::vector<glm::vec<T, n>> &cp) {
    int span = findSpan(u, deg, knots);
    // Insert new control point and knot between span and (span + 1)
    knots.insert(span, u);
    for (int i = 0; i < cp.size(); ++i) {
        cp[i].insert(span, tvecn(0.0));
    }
    // Compute new positions for affected control points
    for (int i = 0; i < cp.size(); ++i) {
        for (int j = span - deg + 1; j <= span; ++j) {
            T a = (u - knots[i]) / (knots[i + deg] - knots[j]);
            cp[i][j] = (1 - a) * cp[i][j - 1] + a * cp[i][j];
        }
    }
}

#endif // REFINE_H
