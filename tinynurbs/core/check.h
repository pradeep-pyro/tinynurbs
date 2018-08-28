#ifndef CHECK_H
#define CHECK_H

#include <vector>
#include <limits>
#include "glm/glm.hpp"
#include "../geometry/curve.h"
#include "../geometry/surface.h"
#include <algorithm>

namespace tinynurbs {

/**
Checks if the relation between degree, number of knots, and
number of control points is valid
@param degree Degree of the NURBS curve
@param num_knots Number of knot values
@param num_ctrl_pts Number of control points
@return Whether the relationship is valid
*/
inline bool isValidRelation(unsigned int degree, size_t num_knots, size_t num_ctrl_pts) {
    return (num_knots - degree - 1) == num_ctrl_pts;
}

template <typename T>
bool isKnotVectorMonotonic(const std::vector<T> &knots) {
    return std::is_sorted(knots.begin(), knots.end());
}

template <typename T>
unsigned int knotMultiplicity(const std::vector<T> &knots, unsigned int index) {
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

template <int dim, typename T>
bool isCurveValid(unsigned int degree, const std::vector<T> &knots, const std::vector<glm::vec<dim, T>> &control_points) {
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

template <int dim, typename T>
bool isCurveValid(unsigned int degree, const std::vector<T> &knots, const std::vector<glm::vec<dim, T>> &control_points,
                  const std::vector<T> &weights) {
    if (!isValid(degree, knots, control_points)) {
        return false;
    }
    if (weights.size() != control_points.size()) {
        return false;
    }
    return true;
}

template <int dim, typename T>
bool isCurveValid(const Curve<dim, T> &crv) {
    return isCurveValid(crv.degree, crv.knots, crv.control_points);
}

template <int dim, typename T>
bool isCurveValid(const RationalCurve<dim, T> &crv) {
    return isCurveValid(crv.degree, crv.knots, crv.control_points, crv.weights);
}

template <int dim, typename T>
bool isSurfaceValid(unsigned int degree_u, unsigned int degree_v, const std::vector<T> &knots_u, const std::vector<T> &knots_v,
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

template <int dim, typename T>
bool isSurfaceValid(unsigned int degree_u, unsigned int degree_v, const std::vector<T> &knots_u,
                    const std::vector<T> &knots_v, const array2<glm::vec<dim, T>> &control_points,
                    const array2<T> &weights) {
    if (!isSurfaceValid(degree_u, degree_v, knots_u, knots_v, control_points)) {
        return false;
    }
    if (control_points.rows() != weights.rows() || control_points.cols() != weights.cols()) {
        return false;
    }
    return true;
}

template <int dim, typename T>
bool isSurfaceValid(const Surface<dim, T> &srf) {
    return isSurfaceValid(srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v, srf.control_points);
}

template <int dim, typename T>
bool isSurfaceValid(const RationalSurface<dim, T> &srf) {
    return isSurfaceValid(srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v, srf.control_points, srf.weights);
}

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

template <int dim, typename T>
bool isCurveClosed(const Curve<dim, T> &crv) {
    return isArray1Closed(crv.degree, crv.control_points) &&
           isKnotVectorClosed(crv.degree, crv.knots);
}

template <int dim, typename T>
bool isCurveClosed(const RationalCurve<dim, T> &crv) {
    return isArray1Closed(crv.degree, crv.control_points) &&
           isArray1Closed(crv.degree, crv.weights) &&
           isKnotVectorClosed(crv.degree, crv.knots);
}

/**
Returns whether the surface is closed along the u-direction
*/
template <int dim, typename T>
bool isSurfaceClosedU(const Surface<dim, T> &srf) {
    return isArray2ClosedU(srf.degree_u, srf.control_points) &&
           isKnotVectorClosed(srf.degree_u, srf.knots_u);
}

/**
Returns whether the surface is closed along the v-direction
*/
template <int dim, typename T>
bool isSurfaceClosedV(const Surface<dim, T> &srf) {
    return isArray2ClosedV(srf.degree_v, srf.control_points) &&
           isKnotVectorClosed(srf.degree_v, srf.knots_v);
}

/**
Returns whether the rational surface is closed along the u-direction
*/
template <int dim, typename T>
bool isSurfaceClosedU(const RationalSurface<dim, T> &srf) {
    return isArray2ClosedU(srf.degree_u, srf.control_points) &&
           isKnotVectorClosed(srf.degree_u, srf.knots_u) &&
           isArray2ClosedU(srf.degree_u, srf.weights);
}

/**
Returns whether the rational surface is closed along the v-direction
*/
template <int dim, typename T>
bool isSurfaceClosedV(const RationalSurface<dim, T> &srf) {
    return isArray2ClosedV(srf.degree_v, srf.control_points) &&
           isKnotVectorClosed(srf.degree_v, srf.knots_v) &&
           isArray2ClosedV(srf.degree_v, srf.weights);
}

} // namespace tinynurbs

#endif // CHECK_H
