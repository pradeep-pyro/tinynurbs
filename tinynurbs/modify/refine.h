#ifndef REFINE_H
#define REFINE_H

#include <vector>
#include "glm/glm.hpp"

namespace tinynurbs {

template <int dim, typename T>
void curveInsertKnot(unsigned int deg, std::vector<T> &knots, std::vector<glm::vec<dim, T>> &cp,
                     T u, unsigned int r=1) {
    int k = findSpan(deg, knots, u);
    int s = 0; // TODO: compute multiplicity
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

} // namespace tinynurbs

#endif // REFINE_H
