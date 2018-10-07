/**
@file
@brief Functions for modifying NURBS curves and surfaces.

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE file.
*/

#ifndef TINYNURBS_MODIFY_H
#define TINYNURBS_MODIFY_H

#include <vector>
#include "glm/glm.hpp"
#include "check.h"
#include "../util/util.h"
#include "../geometry/curve.h"
#include "../geometry/surface.h"

namespace tinynurbs {

/////////////////////////////////////////////////////////////////////

namespace internal {

/**
 * Insert knots in the curve
 * @param deg Degree of the curve
 * @param knots Knot vector of the curve
 * @param cp Control points of the curve
 * @param u Parameter to insert knot(s) at
 * @param r Number of times to insert knot
 * @param[inout] new_knots Updated knot vector
 * @param[inout] new_cp Updated control points
 */
template <int dim, typename T>
void curveKnotInsert(unsigned int deg, const std::vector<T> &knots, const std::vector<glm::vec<dim, T>> &cp, T u,
                     unsigned int r, std::vector<T> &new_knots, std::vector<glm::vec<dim, T>> &new_cp) {
    int k = findSpan(deg, knots, u);
    unsigned int s = knotMultiplicity(knots, k);

    // Insert new knots between span and (span + 1)
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
    for (int j = 1; j <= r; ++j) {
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
}


/**
 * Insert knots in the surface along one direction
 * @param degree Degree of the surface along which to insert knot
 * @param knots Knot vector
 * @param cp 2D array of control points
 * @param knot Knot value to insert
 * @param r Number of times to insert
 * @param along_u Whether inserting along u-direction
 * @param[inout] new_knots Updated knot vector
 * @param[inout] new_cp Updated control points
 */
template <int dim, typename T>
void surfaceKnotInsert(unsigned int degree, const std::vector<T> &knots,
                       const array2<glm::vec<dim, T>> &cp, T knot,
                       unsigned int r, bool along_u,
                       std::vector<T> &new_knots, array2<glm::vec<dim, T>> &new_cp) {
    int span = findSpan(degree, knots, knot);
    unsigned int s = knotMultiplicity(knots, span);

    // Create a new knot vector
    new_knots.resize(knots.size() + r);
    for(int i = 0; i <= span; ++i) {
        new_knots[i] = knots[i];
    }
    for(int i = 1; i <= r; ++i) {
        new_knots[span + i] = knot;
    }
    for(int i = span + 1; i < knots.size(); ++i) {
        new_knots[i + r] = knots[i];
    }
    // Compute alpha
    array2<T> alpha(degree - s, r + 1, T(0));
    for (int j = 1; j <= r; ++j) {
        int L = span - degree + j;
        for (int i = 0; i <= degree - j - s; ++i) {
            alpha(i, j) = (knot - knots[L + i]) / (knots[i + span + 1] - knots[L + i]);
        }
    }

    // Create a temporary container for affected control points per row/column
    std::vector<glm::vec<dim, T>> tmp(degree + 1);

    if (along_u) {
        // Create new control points with additional rows
        new_cp.resize(cp.rows() + r, cp.cols());

        // Update control points
        // Each row is a u-isocurve, each col is a v-isocurve
        for (int col = 0; col < cp.cols(); ++col) {
            // Copy unaffected control points
            for (int i = 0; i <= span - degree; ++i) {
                new_cp(i, col) = cp(i, col);
            }
            for (int i = span - s; i < cp.rows(); ++i) {
                new_cp(i + r, col) = cp(i, col);
            }
            // Copy affected control points to temp array
            for (int i = 0; i < degree - s + 1; ++i) {
                tmp[i] = cp(span - degree + i, col);
            }
            // Insert knot
            for (int j = 1; j <= r; ++j) {
                int L = span - degree + j;
                for (int i = 0; i <= degree - j - s; ++i) {
                    T a = alpha(i, j);
                    tmp[i] = (1 - a) * tmp[i] + a * tmp[i + 1];
                }
                new_cp(L, col) = tmp[0];
                new_cp(span + r - j - s, col) = tmp[degree - j - s];
            }
            int L = span - degree + r;
            for (int i = L + 1; i < span - s; ++i) {
                new_cp(i, col) = tmp[i - L];
            }
        }
    }
    else {
        // Create new control points with additional columns
        new_cp.resize(cp.rows(), cp.cols() + r);

        // Update control points
        // Each row is a u-isocurve, each col is a v-isocurve
        for (int row = 0; row < cp.rows(); ++row) {
            // Copy unaffected control points
            for (int i = 0; i <= span - degree; ++i) {
                new_cp(row, i) = cp(row, i);
            }
            for (int i = span - s; i < cp.cols(); ++i) {
                new_cp(row, i + r) = cp(row, i);
            }
            // Copy affected control points to temp array
            for (int i = 0; i < degree - s + 1; ++i) {
                tmp[i] = cp(row, span - degree + i);
            }
            // Insert knot
            for (int j = 1; j <= r; ++j) {
                int L = span - degree + j;
                for (int i = 0; i <= degree - j - s; ++i) {
                    T a = alpha(i, j);
                    tmp[i] = (1 - a) * tmp[i] + a * tmp[i + 1];
                }
                new_cp(row, L) = tmp[0];
                new_cp(row, span + r - j - s) = tmp[degree - j - s];
            }
            int L = span - degree + r;
            for (int i = L + 1; i < span - s; ++i) {
                new_cp(row, i) = tmp[i - L];
            }
        }
    }
}

/**
 * Unclamp the curve by removing knots and try to retain original shape
 * @param degree Degree of the curve
 * @param knots Knot vector of the curve
 * @param control_points Array of control points of the curve
 * @param start Whether to unclamp at start of curve
 * @param end Whether to unclamp at end of curve
 */
template <int dim, typename T>
void curveUnclamp(unsigned int degree, std::vector<T> &knots,
                  std::vector<glm::vec<dim, T>> control_points,
                  bool start=true, bool end=true) {
    int n = knots.size() - degree - 2;
    int p = degree;
    if (start) {
        for (int i = 0; i < p - 1; ++i)  {
            // Incrementally update knot spacing such the result would be a uniform knot vector when
            // wrapped around
            knots[p - i - 1] = knots[p - i] - (knots[p - i + 1] - knots[n - i]);
            // Update control points to retain shape
            int k = p - 1;
            for (int j = i; j >=0; --j) {
                T alpha = (knots[p] - knots[k]) / (knots[p + j + 1] - knots[k]);
                control_points[j] = (control_points[j] - alpha * control_points[j + 1]) /
                                    (T(1) - alpha);
                --k;
            }
        }
        // Update first knot value
        knots[0] = knots[1] - (knots[n - p + 2] - knots[n - p + 1]);
    }
    if (end) {
        for (int i = 0; i < p - 1; ++i) {
            knots[n + i + 2] = knots[n + i + 1] + (knots[p + i + 1] - knots[p + i]);
            for (int j = i; j >= 0; --j) {
                T alpha = (knots[n + 1] - knots[n + j]) / (knots[n - j + i + 2] - knots[n - j]);
                control_points[j] = (control_points[n - j] - (T(1) - alpha) * control_points[n - j - 1]) / alpha;
            }
        }
        // Update last knot value
        knots[n + p + 1] = knots[n + p] + (knots[2 * p] - knots[2 * p - 1]);
    }
}



} // namespace internal

/////////////////////////////////////////////////////////////////////

/**
 * Insert knots in the curve
 * @param[inout] crv Curve object
 * @param u Parameter to insert knot at
 * @param repeat Number of times to insert
 */
template <int dim, typename T>
void curveKnotInsert(Curve<dim, T> &crv, T u, unsigned int repeat=1) {
    std::vector<T> new_knots;
    std::vector<glm::vec<dim, T>> new_cp;
    internal::curveKnotInsert(crv.degree, crv.knots, crv.control_points, u,
                              repeat, new_knots, new_cp);
    crv.knots = new_knots;
    crv.control_points = new_cp;
}

/**
 * Insert knots in the rational curve
 * @param[inout] crv RationalCurve object
 * @param u Parameter to insert knot at
 * @param repeat Number of times to insert
 */
template <int dim, typename T>
void curveKnotInsert(RationalCurve<dim, T> &crv, T u, unsigned int repeat=1) {
    // Convert to homogenous coordinates
    std::vector<glm::vec<dim + 1, T>> Cw;
    Cw.reserve(crv.control_points.size());
    for (int i = 0; i < crv.control_points.size(); ++i) {
        Cw.push_back(util::cartesianToHomogenous(crv.control_points[i], crv.weights[i]));
    }

    // Perform knot insertion and get new knots and control points
    std::vector<glm::vec<dim + 1, T>> new_Cw;
    std::vector<T> new_knots;
    internal::curveKnotInsert(crv.degree, crv.knots, Cw, u, repeat, new_knots, new_Cw);

    // Convert back to cartesian coordinates
    std::vector<glm::vec<dim, T>> new_cp;
    std::vector<T> new_w;
    new_cp.reserve(new_Cw.size());
    new_w.reserve(new_Cw.size());
    for (int i = 0; i < new_Cw.size(); ++i) {
        new_cp.push_back(util::homogenousToCartesian(new_Cw[i]));
        new_w.push_back(new_Cw[i][dim]);
    }

    // Copy to crv
    crv.knots = new_knots;
    crv.control_points = new_cp;
    crv.weights = new_w;
}

/**
 * Insert knots in the surface along u-direction
 * @param[inout] srf Surface object
 * @param u Knot value to insert
 * @param repeat Number of times to insert
 */
template <int dim, typename T>
void surfaceKnotInsertU(Surface<dim, T> &srf, T u, unsigned int repeat=1) {
    // New knots and new control points after knot insertion
    std::vector<T> new_knots_u;
    array2<glm::vec<dim, T>> new_cp;
    surfaceKnotInsert(srf.degree_u, srf.knots_u, srf.control_points, u, repeat, true,
                      new_knots_u, new_cp);
    // Copy to given surface
    srf.knots_u = new_knots_u;
    srf.control_points = new_cp;
}

/**
 * Insert knots in the rational surface along u-direction
 * @param[inout] srf RationalSurface object
 * @param u Knot value to insert
 * @param repeat Number of times to insert
 */
template <int dim, typename T>
void surfaceKnotInsertU(RationalSurface<dim, T> &srf, T u, unsigned int repeat=1) {
    // Original control points in homogenous coordinates
    array2<glm::vec<dim + 1, T>> Cw(srf.control_points.rows(), srf.control_points.cols());
    for (int i = 0; i < srf.control_points.rows(); ++i) {
        for (int j = 0; j < srf.control_points.cols(); ++j) {
            Cw(i, j) = util::cartesianToHomogenous(srf.control_points(i, j), srf.weights(i, j));
        }
    }

    // New knots and new homogenous control points after knot insertion
    std::vector<T> new_knots_u;
    array2<glm::vec<dim + 1, T>> new_Cw;
    surfaceKnotInsert(srf.degree_u, srf.knots_u, Cw, u, repeat, true,
                      new_knots_u, new_Cw);

    // Convert back to cartesian coordinates
    array2<glm::vec<dim, T>> new_cp(new_Cw.rows(), new_Cw.cols());
    array2<T> new_w(new_Cw.rows(), new_Cw.cols());
    for (int i = 0; i < new_Cw.rows(); ++i) {
        for (int j = 0; j < new_Cw.cols(); ++j) {
            new_cp(i, j) = util::homogenousToCartesian(new_Cw(i, j));
            new_w(i, j) = new_Cw(i, j)[dim];
        }
    }
    srf.knots_u = new_knots_u;
    srf.control_points = new_cp;
    srf.weights = new_w;
}

/**
 * Insert knots in the surface along v-direction
 * @param[inout] srf Surface object
 * @param v Knot value to insert
 * @param repeat Number of times to insert
 */
template <int dim, typename T>
void surfaceKnotInsertV(Surface<dim, T> &srf, T v, unsigned int repeat=1) {
    // New knots and new control points after knot insertion
    std::vector<T> new_knots_v;
    array2<glm::vec<dim, T>> new_cp;
    surfaceKnotInsert(srf.degree_u, srf.knots_u, srf.control_points, v, repeat, false,
                      new_knots_v, new_cp);
    // Copy to given surface
    srf.knots_v = new_knots_v;
    srf.control_points = new_cp;
}

/**
 * Insert knots in the rational surface along v-direction
 * @param[inout] srf RationalSurface object
 * @param v Knot value to insert
 * @param repeat Number of times to insert
 */
template <int dim, typename T>
void surfaceKnotInsertV(RationalSurface<dim, T> &srf, T v, unsigned int repeat=1) {
    // Original control points in homogenous coordinates
    array2<glm::vec<dim + 1, T>> Cw(srf.control_points.rows(), srf.control_points.cols());
    for (int i = 0; i < srf.control_points.rows(); ++i) {
        for (int j = 0; j < srf.control_points.cols(); ++j) {
            Cw(i, j) = util::cartesianToHomogenous(srf.control_points(i, j), srf.weights(i, j));
        }
    }

    // New knots and new homogenous control points after knot insertion
    std::vector<T> new_knots_v;
    array2<glm::vec<dim + 1, T>> new_Cw;
    surfaceKnotInsert(srf.degree_u, srf.knots_u, Cw, v, repeat, false,
                      new_knots_v, new_Cw);

    // Convert back to cartesian coordinates
    array2<glm::vec<dim, T>> new_cp(new_Cw.rows(), new_Cw.cols());
    array2<T> new_w(new_Cw.rows(), new_Cw.cols());
    for (int i = 0; i < new_Cw.rows(); ++i) {
        for (int j = 0; j < new_Cw.cols(); ++j) {
            new_cp(i, j) = util::homogenousToCartesian(new_Cw(i, j));
            new_w(i, j) = new_Cw(i, j)[dim];
        }
    }
    srf.knots_v = new_knots_v;
    srf.control_points = new_cp;
    srf.weights = new_w;
}

/**
 * Unclamp the curve by removing knots and trying to retain original shape
 * @param crv Curve object
 * @param start Whether to unclamp at start of curve
 * @param end Whether to unclamp at end of curve
 */
template <int dim, typename T>
void curveUnclamp(Curve<dim, T> &crv, bool start=true, bool end=true) {
    curveUnclamp(crv.degree, crv.knots, crv.control_points, start, end);
}

/**
 * Unclamp the rational curve by removing knots and trying to retain original shape
 * @param crv RationalCurve object
 * @param start Whether to unclamp at start of curve
 * @param end Whether to unclamp at end of curve
 */
template <int dim, typename T>
void curveUnclamp(RationalCurve<dim, T> &crv, bool start=true, bool end=true) {
    std::vector<glm::vec<dim + 1, T>> Cw;
    Cw.reserve(crv.control_points.size());
    for (int i = 0; i < crv.control_points.size(); ++i) {
        Cw.push_back(util::cartesianToHomogenous(crv.control_points[i], crv.weights[i]));
    }
    curveUnclamp(crv.degree, crv.knots, Cw, start, end);
    for (int i = 0; i < crv.control_points.size(); ++i) {
        crv.control_points[i] = util::homogenousToCartesian(Cw[i]);
        crv.weights[i] = Cw[dim];
    }
}

/**
 * Clamp the curve by adding knots while retaining original shape
 * @param crv Curve object
 * @param start Whether to clamp at start of curve
 * @param end Whether to clamp at end of curve
 */
template <int dim, typename T>
void curveClamp(Curve<dim, T> &crv, bool start=true, bool end=true) {
    if (start) {
        T u = crv.knots[crv.degree];
        curveKnotInsert(crv, u, crv.degree - 1);
    }
    if (end) {
        T u = crv.knots[crv.knots.size() - crv.degree - 1];
        curveKnotInsert(crv, u, crv.degree - 1);
    }
}

/**
 * Clamp the rational curve by adding knots while retaining original shape
 * @param crv RationalCurve object
 * @param start Whether to clamp at start of curve
 * @param end Whether to clamp at end of curve
 */
template <int dim, typename T>
void curveClamp(RationalCurve<dim, T> &crv, bool start=true, bool end=true) {
    if (start) {
        T u = crv.knots[crv.degree];
        curveKnotInsert(crv, u, crv.degree - 1);
    }
    if (end) {
        T u = crv.knots[crv.knots.size() - crv.degree - 1];
        curveKnotInsert(crv, u, crv.degree - 1);
    }
}

} // namespace tinynurbs

#endif // TINYNURBS_MODIFY_H
