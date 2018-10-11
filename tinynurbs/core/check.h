/**
@file
@brief Functionality for checking validity and properties of NURBS curves and
surfaces

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE file.
*/

#ifndef TINYNURBS_CHECK_H
#define TINYNURBS_CHECK_H

#include <vector>
#include <limits>
#include "glm/glm.hpp"
#include "curve.h"
#include "surface.h"
#include <algorithm>

namespace tinynurbs {

/////////////////////////////////////////////////////////////////////

namespace internal {

/**
Checks if the relation between degree, number of knots, and
number of control points is valid
@param degree Degree
@param num_knots Number of knot values
@param num_ctrl_pts Number of control points
@return Whether the relationship is valid
*/
inline bool isValidRelation(unsigned int degree, size_t num_knots,
                            size_t num_ctrl_pts) {
    return (num_knots - degree - 1) == num_ctrl_pts;
}

/**
 * isKnotVectorMonotonic returns whether the knots are in ascending order
 * @tparam Type of knot values
 * @param knots Knot vector
 * @return Whether monotonic
 */
template <typename T>
bool isKnotVectorMonotonic(const std::vector<T> &knots) {
    return std::is_sorted(knots.begin(), knots.end());
}

/**
 * Returns whether the curve is valid
 * @tparam dim Dimension of the curve (2 or 3 in general)
 * @tparam T Type of control point coordinates, knot values
 * @param degree Degree of curve
 * @param knots Knot vector of curve
 * @param control_points Control points of curve
 * @return Whether valid
 */
template <int dim, typename T>
bool curveIsValid(unsigned int degree, const std::vector<T> &knots,
                  const std::vector<glm::vec<dim, T>> &control_points) {
    if (degree < 1 || degree > 9) {
        return false;
    }
    if (!isValidRelation(degree, knots.size(), control_points.size())) {
        return false;
    }
    if (!isKnotVectorMonotonic(knots)) {
        return false;
    }
    return true;
}

/**
 * Returns whether the curve is valid
 * @tparam dim Dimension of the curve (2 or 3 in general)
 * @tparam T Type of control point coordinates, knot values and weights
 * @param degree Degree of curve
 * @param knots Knot vector of curve
 * @param control_points Control points of curve
 * @return Whether valid
 */
template <int dim, typename T>
bool curveIsValid(unsigned int degree, const std::vector<T> &knots,
                  const std::vector<glm::vec<dim, T>> &control_points,
                  const std::vector<T> &weights) {
    if (!isValidRelation(degree, knots.size(), control_points.size())) {
        return false;
    }
    if (weights.size() != control_points.size()) {
        return false;
    }
    return true;
}

/**
 * Returns whether the surface is valid
 * @tparam dim Dimension of the surface (3 in general)
 * @tparam T Type of control point coordinates, knot values
 * @param degree_u Degree of surface along u-direction
 * @param degree_v Degree of surface along v-direction
 * @param knots_u Knot vector of surface along u-direction
 * @param knots_v Knot vector of surface along v-direction
 * @param control_points Control points grid of surface
 * @return Whether valid
 */
template <int dim, typename T>
bool surfaceIsValid(unsigned int degree_u, unsigned int degree_v,
                    const std::vector<T> &knots_u, const std::vector<T> &knots_v,
                    const array2<glm::vec<dim, T>> &control_points) {
    if (degree_u < 1 || degree_u > 9 || degree_v < 1 || degree_v > 9) {
        return false;
    }
    if (!isValidRelation(degree_u, knots_u.size(), control_points.rows()) ||
            !isValidRelation(degree_v, knots_v.size(), control_points.cols())) {
        return false;
    }
    if (!isKnotVectorMonotonic(knots_u) || !isKnotVectorMonotonic(knots_v)) {
        return false;
    }
    return true;
}

/**
 * Returns whether the rational surface is valid
 * @tparam dim Dimension of the surface (3 in general)
 * @tparam T Type of control point coordinates, knot values
 * @param degree_u Degree of surface along u-direction
 * @param degree_v Degree of surface along v-direction
 * @param knots_u Knot vector of surface along u-direction
 * @param knots_v Knot vector of surface along v-direction
 * @param control_points Control points grid of surface
 * @param weights Weights corresponding to control point grid of surface
 * @return Whether valid
 */
template <int dim, typename T>
bool surfaceIsValid(unsigned int degree_u, unsigned int degree_v,
                    const std::vector<T> &knots_u, const std::vector<T> &knots_v,
                    const array2<glm::vec<dim, T>> &control_points,
                    const array2<T> &weights) {
    if (!surfaceIsValid(degree_u, degree_v, knots_u, knots_v, control_points)) {
        return false;
    }
    if (control_points.rows() != weights.rows() ||
            control_points.cols() != weights.cols()) {
        return false;
    }
    return true;
}


/**
 * Returns whether the given knot vector is closed by checking the
 * periodicity of knot vectors near the start and end
 * @param degree Degree of curve/surface
 * @param knots Knot vector of curve/surface
 * @return Whether knot vector is closed
 */
template <typename T>
bool isKnotVectorClosed(unsigned int degree, const std::vector<T> &knots) {
    for (int i = 0; i < degree - 1; ++i) {
        int j = knots.size() - degree + i;
        if (std::abs((knots[i + 1] - knots[i]) - (knots[j + 1] - knots[j])) >
                std::numeric_limits<T>::epsilon()) {
            return false;
        }
    }
    return true;
}

/**
 * Returns whether the given knot vector is closed by checking the
 * periodicity of knot vectors near the start and end
 * @param degree Degree of curve/surface
 * @param vec Array of any control points/weights
 * @return Whether knot vector is closed
 */
template <typename T>
bool isArray1Closed(unsigned int degree, const std::vector<T> &vec) {
    for (int i = 0; i < degree; ++i) {
        int j = vec.size() - degree + i;
        if (glm::length(vec[i] - vec[j]) > 1e-5) {
            return false;
        }
    }
    return true;
}

/**
 * Returns whether the 2D array is closed along the u-direction
 * i.e., along rows.
 * @param degree_u Degree along u-direction
 * @param arr 2D array of control points / weights
 * @return Whether closed along u-direction
 */
template <typename T>
bool isArray2ClosedU(unsigned int degree_u, const array2<T> &arr) {
    for (int i = 0; i < degree_u; ++i) {
        for (int j = 0; j < arr.cols(); ++j) {
            int k = arr.cols() - degree_u + i;
            if (glm::length(arr(i, j) - arr(k, j)) > 1e-5) {
                return false;
            }
        }
    }
    return true;
}

/**
 * Returns whether the 2D array is closed along the v-direction
 * i.e., along columns.
 * @param degree_v Degree along v-direction
 * @param arr 2D array of control points / weights
 * @return Whether closed along v-direction
 */
template <typename T>
bool isArray2ClosedV(unsigned int degree_v, const array2<T> &arr) {
    for (int i = 0; i < arr.rows(); ++i) {
        for (int j = 0; j < degree_v; j++) {
            int k = arr.rows() - degree_v + i;
            if (glm::length(arr(i, j) - arr(i, k)) > 1e-5) {
                return false;
            }
        }
    }
    return true;
}

} // namespace internal

/////////////////////////////////////////////////////////////////////

/**
 * Returns the mulitplicity of the knot at index
 * @tparam Type of knot values
 * @param knots Knot vector
 * @param index Index of knot of interest
 * @return Multiplicity (>= 1)
 */
template <typename T>
unsigned int knotMultiplicity(const std::vector<T> &knots,
                              unsigned int index) {
    T u = knots[index];
    T eps = std::numeric_limits<T>::epsilon();
    unsigned int mult = 0;
    for (unsigned int i = index; i < knots.size(); ++i) {
        if (std::abs(u - knots[index + 1]) < eps) {
            ++mult;
        }
    }
    return mult;
}

/**
 * Returns whether the curve is valid
 * @tparam dim Dimension of the curve (2 or 3 in general)
 * @tparam T Type of control point coordinates, knot values
 * @param crv Curve object
 * @return Whether valid
 */
template <int dim, typename T>
bool curveIsValid(const Curve<dim, T> &crv) {
    return internal::curveIsValid(crv.degree, crv.knots, crv.control_points);
}

/**
 * Returns whether the curve is valid
 * @tparam dim Dimension of the curve (2 or 3 in general)
 * @tparam T Type of control point coordinates, knot values
 * @param crv RationalCurve object
 * @return Whether valid
 */
template <int dim, typename T>
bool curveIsValid(const RationalCurve<dim, T> &crv) {
    return internal::curveIsValid(crv.degree, crv.knots, crv.control_points,
                                  crv.weights);
}

/**
 * Returns whether the surface is valid
 * @tparam dim Dimension of the surface (3 in general)
 * @tparam T Type of control point coordinates, knot values
 * @param srf Surface object
 * @return Whether valid
 */
template <int dim, typename T>
bool surfaceIsValid(const Surface<dim, T> &srf) {
    return internal::surfaceIsValid(srf.degree_u, srf.degree_v, srf.knots_u,
                                    srf.knots_v, srf.control_points);
}

/**
 * Returns whether the rational surface is valid
 * @tparam dim Dimension of the surface (3 in general)
 * @tparam T Type of control point coordinates, knot values
 * @param srf RationalSurface object
 * @return Whether valid
 */
template <int dim, typename T>
bool surfaceIsValid(const RationalSurface<dim, T> &srf) {
    return internal::surfaceIsValid(srf.degree_u, srf.degree_v, srf.knots_u,
                                    srf.knots_v, srf.control_points, srf.weights);
}

/**
 * Checks whether the curve is closed
 * @param crv Curve object
 * @return  Whether closed
 */
template <int dim, typename T>
bool curveIsClosed(const Curve<dim, T> &crv) {
    return internal::isArray1Closed(crv.degree, crv.control_points) &&
           internal::isKnotVectorClosed(crv.degree, crv.knots);
}

/**
 * Checks whether the rational curve is closed
 * @param crv RationalCurve object
 * @return  Whether closed
 */
template <int dim, typename T>
bool curveIsClosed(const RationalCurve<dim, T> &crv) {
    return internal::isArray1Closed(crv.degree, crv.control_points) &&
           internal::isArray1Closed(crv.degree, crv.weights) &&
           internal::isKnotVectorClosed(crv.degree, crv.knots);
}

/**
 * Checks whether the surface is closed along u-direction
 * @param srf Surface object
 * @return  Whether closed along u-direction
*/
template <int dim, typename T>
bool surfaceIsClosedU(const Surface<dim, T> &srf) {
    return internal::isArray2ClosedU(srf.degree_u, srf.control_points) &&
           internal::isKnotVectorClosed(srf.degree_u, srf.knots_u);
}

/**
 * Checks whether the surface is closed along v-direction
 * @param srf Surface object
 * @return  Whether closed along v-direction
*/
template <int dim, typename T>
bool surfaceIsClosedV(const Surface<dim, T> &srf) {
    return internal::isArray2ClosedV(srf.degree_v, srf.control_points) &&
           internal::isKnotVectorClosed(srf.degree_v, srf.knots_v);
}

/**
 * Checks whether the rational surface is closed along u-direction
 * @param srf RationalSurface object
 * @return  Whether closed along u-direction
*/
template <int dim, typename T>
bool surfaceIsClosedU(const RationalSurface<dim, T> &srf) {
    return internal::isArray2ClosedU(srf.degree_u, srf.control_points) &&
           internal::isKnotVectorClosed(srf.degree_u, srf.knots_u) &&
           internal::isArray2ClosedU(srf.degree_u, srf.weights);
}

/**
 * Checks whether the rational surface is closed along v-direction
 * @param srf RationalSurface object
 * @return  Whether closed along v-direction
*/
template <int dim, typename T>
bool surfaceIsClosedV(const RationalSurface<dim, T> &srf) {
    return internal::isArray2ClosedV(srf.degree_v, srf.control_points) &&
           internal::isKnotVectorClosed(srf.degree_v, srf.knots_v) &&
           internal::isArray2ClosedV(srf.degree_v, srf.weights);
}

} // namespace tinynurbs

#endif // TINYNURBS_CHECK_H
