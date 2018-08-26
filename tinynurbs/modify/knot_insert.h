#ifndef REFINE_H
#define REFINE_H

#include <vector>
#include "glm/glm.hpp"
#include "../core/check.h"
#include "../util/util.h"
#include "../geometry/curve.h"

namespace tinynurbs {

template <int dim, typename T>
void curveKnotInsert(unsigned int deg, std::vector<T> &knots, std::vector<glm::vec<dim, T>> &cp,
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
void curveKnotInsert(Curve<dim, T> &crv, T u, unsigned int repeat=1) {
    curveKnotInsert(crv.degree, crv.knots, crv.control_points, u, repeat);
}

template <int dim, typename T>
void curveKnotInsert(RationalCurve<dim, T> &crv, T u, unsigned int repeat=1) {
    std::vector<glm::vec<dim + 1, T>> Cw;
    Cw.reserve(crv.control_points.size());
    for (int i = 0; i < crv.control_points.size(); ++i) {
        Cw.push_back(util::cartesianToHomogenous(crv.control_points[i], crv.weights[i]));
    }
    curveKnotInsert(crv.degree, crv.knots, Cw, u, repeat);
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

template <int dim, typename T>
void surfaceKnotInsertU(unsigned int deg_u, unsigned int deg_v,
                        std::vector<T> &knots_u, std::vector<T> &knots_v,
                        array2<glm::vec<dim, T>> &cp,
                        T u, unsigned int r=1) {
    int span = findSpan(deg_u, knots_u, u);
    unsigned int s = knotMultiplicity(knots_u, span);

    // Create a new knot vector
    std::vector<T> new_knots_u;
    new_knots_u.resize(knots_u.size() + r);
    for(int i = 0; i < span + 1; ++i) {
        new_knots_u[i] = knots_u[i];
    }
    for(int i = 1; i < r + 1; ++i) {
        new_knots_u[span + i] = u;
    }
    for(int i = span + 1; i < knots_u.size(); ++i) {
        new_knots_u[i + r] = knots_u[i];
    }
    // Create new control points
    array2<glm::vec<dim, T>> new_cp(cp.rows() + r, cp.cols());
    // Create a temporary container for affected control points per row
    std::vector<glm::vec<dim, T>> tmp_cp(deg_u + 1);

    // Compute alpha
    array2<T> alpha(r, deg_u - s, T(0));
    for (int j = 1; j <= r; ++j) {
        int L = span - deg_u + j;
        for (int i = 0; i <= deg_u - j - s; ++i) {
            alpha(i, j) = (u - knots_u[L + i]) / (knots_u[i + span + 1] - knots_u[L + i]);
        }
    }
    // Update control points
    for (int row = 0; row < new_cp.rows(); ++row) {
        // Copy unaffected control points
        for (int col = 0; col < span - deg_u + 1; ++col) {
            new_cp(row, col) = cp(row, col);
        }
        for (int col = span; col < new_cp.cols(); ++col) {
            new_cp(row, col + r) = cp(row, col);
        }
        // Copy affected control points to temp array
        for (int i = 0; i < deg_u - s + 1; ++i) {
            tmp_cp[i] = cp(row, span - deg_u + i);
        }
        // Insert knot
        for (int j = 1; j <= r; ++j) {
            int L = span - deg_u + j;
            for (int i = 0; i <= deg_u - j - s; ++i) {
                T a = alpha(i, j);
                tmp_cp[i] = a * tmp_cp[i + 1] + (T(1) - a) * tmp_cp[i];
            }
            new_cp(row, L) = tmp_cp[0];
            new_cp(row, span + r - j) = tmp_cp[deg_u - j];
        }
        int L = span - deg_u + r;
        for (int i = L + 1; i < span; ++i) {
            new_cp(row, i) = tmp_cp[i - L];
        }
    }
    cp = new_cp;
    knots_u = new_knots_u;
}

} // namespace tinynurbs

#endif // REFINE_H
