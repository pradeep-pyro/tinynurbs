#ifndef CHECK_H
#define CHECK_H

#include <vector>
#include <limits>
#include "glm/glm.hpp"
#include "../geometry/curve.h"
#include <algorithm>

namespace tinynurbs {

template <typename T>
bool isKnotVectorMonotonic(const std::vector<T> &knots) {
    return std::is_sorted(knots.begin(), knots.end());
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
bool isRationalCurveValid(unsigned int degree, const std::vector<T> &knots, const std::vector<glm::vec<dim, T>> &control_points,
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
bool isRationalCurveValid(const RationalCurve<dim, T> &crv) {
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
bool isRationalSurfaceValid(unsigned int degree_u, unsigned int degree_v, const std::vector<T> &knots_u,
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
bool isRationalSurfaceValid(const RationalSurface<dim, T> &srf) {
    return isRationalSurfaceValid(srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v, srf.control_points, srf.weights);
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

template <int dim, typename T>
bool isControlPointsClosed(unsigned int degree, const std::vector<glm::vec<dim, T>> &control_points) {
    for (int i = 0; i < degree; ++i) {
        int j = control_points.size() - degree + i;
        if (glm::length(control_points[i] - control_points[j]) >
                std::numeric_limits<T>::epsilon()) {
            return false;
        }
    }
    return true;
}

template <int dim, typename T>
bool isWeightsClosed(unsigned int degree, const std::vector<T> &weights) {
    for (int i = 0; i < degree; ++i) {
        int j = weights.size() - degree + i;
        if (std::abs(weights[i] - weights[j]) > std::numeric_limits<T>::epsilon()) {
            return false;
        }
    }
    return true;
}

template <int dim, typename T>
bool isControlPointsClosedU(unsigned int degree_u, const array2<glm::vec<dim, T>> &control_points) {
    for (int i = 0; i < degree_u; ++i) {
        for (int j = 0; j < control_points.cols(); ++j) {
            int k = control_points.cols() - degree_u + i;
            if (glm::length(control_points(i, j) - control_points(k, j)) >
                    std::numeric_limits<T>::epsilon()) {
                return false;
            }
        }
    }
    return true;
}

template <int dim, typename T>
bool isControlPointsClosedV(unsigned int degree_v, const array2<glm::vec<dim, T>> &control_points) {
    for (int i = 0; i < control_points.rows(); ++i) {
        for (int j = 0; j < degree_v; j++) {
            int k = control_points.rows() - degree_v + i;
            if (glm::length(control_points(i, j) - control_points(i, k)) >
                    std::numeric_limits<T>::epsilon()) {
                return false;
            }
        }
    }
    return true;
}

template <typename T>
bool isWeightsClosedU(unsigned int degree_u, const array2<T> &weights) {
    for (int i = 0; i < degree_u; ++i) {
        for (int j = 0; j < weights.rows(); ++j) {
            int k = weights.rows() - degree_u + i;
            if (std::abs(weights(i, j) - weights(k, j)) > std::numeric_limits<T>::epsilon()) {
                return false;
            }
        }
    }
    return true;
}

template <typename T>
bool isWeightsClosedV(unsigned int degree_v, const array2<T> &weights) {
    for (int i = 0; i < weights.rows(); ++i) {
        for (int j = 0; j < degree_v; j++) {
            int k = weights.rows() - degree_v + i;
            if (glm::length(weights(i, j) - weights(i, k)) >
                    std::numeric_limits<T>::epsilon()) {
                return false;
            }
        }
    }
    return true;
}


template <int dim, typename T>
bool isCurveClosed(const Curve<dim, T> &crv) {
    return isControlPointsClosed(crv.degree, crv.control_points) &&
           isKnotVectorClosed(crv.degree, crv.knots);
}

template <int dim, typename T>
bool isRationalCurveClosed(const RationalCurve<dim, T> &crv) {
    return isControlPointsClosed(crv.degree, crv.control_points) &&
           isWeightsClosed(crv.degree, crv.weights) &&
           isKnotVectorClosed(crv.degree, crv.knots);
}

/**
Returns whether the surface is closed along the u-direction
*/
template <int dim, typename T>
bool isSurfaceClosedU(const Surface<dim, T> &srf) {
    return isControlPointsClosedU(srf.degree_u, srf.control_points) &&
           isKnotVectorClosed(srf.degree_u, srf.knots_u);
}

/**
Returns whether the surface is closed along the v-direction
*/
template <int dim, typename T>
bool isSurfaceClosedV(const Surface<dim, T> &srf) {
    return isControlPointsClosedV(srf.degree_v, srf.control_points) &&
           isKnotVectorClosed(srf.degree_v, srf.knots_v);
}

/**
Returns whether the rational surface is closed along the u-direction
*/
template <int dim, typename T>
bool isRationalSurfaceClosedU(const RationalSurface<dim, T> &srf) {
    return isControlPointsClosedU(srf.degree_u, srf.control_points) &&
           isKnotVectorClosed(srf.degree_u, srf.knots_u) &&
           isWeightsClosedU(srf.degree_u, srf.weights);
}

/**
Returns whether the rational surface is closed along the v-direction
*/
template <int dim, typename T>
bool isRationalSurfaceClosedV(const RationalSurface<dim, T> &srf) {
    return isControlPointsClosedV(srf.degree_v, srf.control_points) &&
           isKnotVectorClosed(srf.degree_v, srf.knots_v) &&
           isWeightsClosedV(srf.degree_v, srf.weights);
}

} // namespace tinynurbs

#endif // CHECK_H
