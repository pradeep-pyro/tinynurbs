/*
@file nurbstk/evaluate.h
@author Pradeep Kumar Jayaraman <pradeep.pyro@gmail.com>

Core functionality for evaluating points and derivatives on NURBS curves and
surfaces

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE.txt file.
*/

#pragma once

#include <vector>

#include "glm/glm.hpp"

#include "basis.h"
#include "array2.h"
#include "util.h"
#include "curve.h"
#include "surface.h"

namespace nurbstk {

/**
Evaluate point on a nonrational NURBS curve
@param[in] degree Degree of the given curve.
@param[in] knots Knot vector of the curve.
@param[in] control_points Control points of the curve.
@param[in] u Parameter to evaluate the curve at.
@param[in, out] point Resulting point on the curve at parameter u.
*/
template <int dim, typename T>
void curvePoint(unsigned int degree, const std::vector<T> &knots,
                const std::vector<glm::vec<dim, T>> &control_points,
                T u, glm::vec<dim, T> &point) {
    // Initialize result to 0s
    for (int i = 0; i < dim; i++) {
        point[i] = static_cast<T>(0.0);
    }

    // Find span and corresponding non-zero basis functions
    int span = findSpan(degree, knots, u);
    std::vector<T> N;
    bsplineBasis(degree, span, knots, u, N);

    // Compute point
    for (int j = 0; j <= degree; j++) {
        point += static_cast<T>(N[j]) * control_points[span - degree + j];
    }
}

/**
Evaluate point on a nonrational NURBS curve
@param[in] crv Curve object
@param[in] u Parameter to evaluate the curve at.
@param[in, out] point Resulting point on the curve at parameter u.
*/
template <int dim, typename T>
void curvePoint(const Curve<dim, T> &crv, T u, glm::vec<dim, T> &point) {
    curvePoint(crv.degree, crv.knots, crv.control_points, u, point);
}

/**
Evaluate point on a rational NURBS curve
@param[in] degree Degree of the curve
@param[in] knots Knot vector of the curve.
@param[in] control_points Control points of the curve.
@param[in] weights Weights corresponding to each control point.
@param[in] u Parameter to evaluate the curve at.
@param[in, out] point Resulting point on the curve.
*/
template <int dim, typename T>
void rationalCurvePoint(unsigned int degree, const std::vector<T> &knots,
                        const std::vector<glm::vec<dim, T>> &control_points,
                        const std::vector<T> &weights,
                        T u, glm::vec<dim, T> &point) {

    typedef glm::vec<dim + 1, T> tvecnp1;

    // Compute homogenous coordinates of control points
    std::vector<tvecnp1> Cw;
    Cw.reserve(control_points.size());
    for (int i = 0; i < control_points.size(); i++) {
        Cw.push_back(tvecnp1(
                         util::cartesianToHomogenous(control_points[i], weights[i])
                     ));
    }

    // Compute point using homogenous coordinates
    tvecnp1 pointw;
    curvePoint(degree, knots, Cw, u, pointw);

    // Convert back to cartesian coordinates
    point = util::homogenousToCartesian(pointw);
}

/**
Evaluate point on a rational NURBS curve
@param[in] crv RationalCurve object
@param[in] u Parameter to evaluate the curve at.
@param[in, out] point Resulting point on the curve.
*/
template <int dim, typename T>
void rationalCurvePoint(const RationalCurve<dim, T> &crv, T u, glm::vec<dim, T> &point) {
    rationalCurvePoint(crv.degree, crv.knots, crv.control_points, crv.weights, u, point);
}

/**
Evaluate derivatives of a non-rational NURBS curve
@param[in] degree Degree of the curve
@param[in] knots Knot vector of the curve.
@param[in] control_points Control points of the curve.
@param[in] num_ders Number of times to derivate.
@param[in] u Parameter to evaluate the derivatives at.
@param[in, out] curve_ders Derivatives of the curve at u.
E.g. curve_ders[n] is the nth derivative at u, where 0 <= n <= num_ders.
*/
template <int dim, typename T>
void curveDerivatives(unsigned int degree, const std::vector<T> &knots,
                      const std::vector<glm::vec<dim, T>> &control_points,
                      int num_ders, T u, std::vector<glm::vec<dim, T>> &curve_ders) {

    typedef glm::vec<dim, T> tvecn;
    using std::vector;

    curve_ders.clear();
    curve_ders.resize(num_ders + 1);

    // Assign higher order derivatives to zero
    for (int k = degree + 1; k <= num_ders; k++) {
        curve_ders[k] = tvecn(0.0);
    }

    // Find the span and corresponding non-zero basis functions & derivatives
    int span = findSpan(degree, knots, u);
    array2<T> ders;
    bsplineDerBasis<T>(degree, span, knots, u, num_ders, ders);

    // Compute first num_ders derivatives
    int du = num_ders < degree ? num_ders : degree;
    for (int k = 0; k <= du; k++) {
        curve_ders[k] = tvecn(0.0);
        for (int j = 0; j <= degree; j++) {
            curve_ders[k] += static_cast<T>(ders(k, j)) *
                            control_points[span - degree + j];
        }
    }
}

/**
Evaluate derivatives of a non-rational NURBS curve
@param[in] crv Curve object
@param[in] num_ders Number of times to derivate.
@param[in] u Parameter to evaluate the derivatives at.
@param[in, out] curve_ders Derivatives of the curve at u.
E.g. curve_ders[n] is the nth derivative at u, where 0 <= n <= num_ders.
*/
template <int dim, typename T>
void curveDerivatives(const Curve<dim, T> &crv, int num_ders, T u,
                      std::vector<glm::vec<dim, T>> &curve_ders) {
    curveDerivatives(crv.degree, crv.knots, crv.control_points, num_ders, u, curve_ders);
}

/**
Evaluate derivatives of a rational NURBS curve
@param[in] u Parameter to evaluate the derivatives at.
@param[in] knots Knot vector of the curve.
@param[in] control_points Control points of the curve.
@param[in] weights Weights corresponding to each control point.
@param[in] num_ders Number of times to differentiate.
@param[in, out] curve_ders Derivatives of the curve at u.
E.g. curve_ders[n] is the nth derivative at u, where n is between 0 and num_ders-1.
*/
template <int dim, typename T>
void rationalCurveDerivatives(unsigned int degree, const std::vector<T> &knots,
                              const std::vector<glm::vec<dim, T>> &control_points,
                              const std::vector<T> weights, int num_ders,
                              T u, std::vector<glm::vec<dim, T>> &curve_ders) {

    typedef glm::vec<dim, T> tvecn;
    typedef glm::vec<dim + 1, T> tvecnp1;
    using std::vector;

    curve_ders.clear();
    curve_ders.resize(num_ders + 1);

    // Compute homogenous coordinates of control points
    vector<tvecnp1> Cw;
    Cw.reserve(control_points.size());
    for (int i = 0; i < control_points.size(); i++) {
        Cw.push_back(tvecnp1(
                         util::cartesianToHomogenous(control_points[i], weights[i])
                     ));
    }

    // Derivatives of Cw
    vector<tvecnp1> Cwders;
    curveDerivatives(u, degree, knots, Cw, num_ders, Cwders);

    // Split Cwders into coordinates and weights
    vector<tvecn> Aders;
    vector<T> wders;
    for (const auto &val : Cwders) {
        tvecn Aderi;
        for (int i = 0; i < dim - 1; i++) {
            Aderi[i] = val[i];
        }
        Aders.push_back(Aderi);
        wders.push_back(val[dim - 1]);
    }

    // Compute rational derivatives
    for (int k = 0; k <= num_ders; k++) {
        tvecn v = Aders[k];
        for (int i = 1; i <= k; i++) {
            v -= static_cast<T>(util::binomial(k, i)) * wders[i] * curve_ders[k - i];
        }
        curve_ders[k] = v / wders[0];
    }
}

/**
Evaluate derivatives of a rational NURBS curve
@param[in] crv RationalCurve object
@param[in] num_ders Number of times to differentiate.
@param[in, out] curve_ders Derivatives of the curve at u.
E.g. curve_ders[n] is the nth derivative at u, where n is between 0 and num_ders-1.
*/
template <int dim, typename T>
void rationalCurveDerivatives(const RationalCurve<dim, T> &crv, int num_ders,
                              T u, std::vector<glm::vec<dim, T>> &curve_ders) {
    rationalCurveDerivatives(crv.degree, crv.knots, crv.control_points, crv.weights,
                             num_ders, u, curve_ders);
}

/**
Evaluate the tangent of a B-spline curve
@param[in] crv Curve object
@param[in, out] tgt Unit tangent of the curve at u.
*/
template <int dim, typename T>
void curveTangent(const Curve<dim, T> &crv, T u, glm::vec<dim, T> &tgt) {
    std::vector<glm::vec<dim, T>> &ders;
    curveDerivatives(crv, 1, u, ders);
    tgt = glm::normalize(ders[1]);
}

/**
Evaluate the tangent of a rational B-spline curve
@param[in] crv RationalCurve object
@param[in, out] tgt Unit tangent of the curve at u.
*/
template <int dim, typename T>
void rationalCurveTangent(const RationalCurve<dim, T> &crv, T u, glm::vec<dim, T> &tgt) {
    std::vector<glm::vec<dim, T>> &ders;
    rationalCurveDerivatives(crv, 1, u, ders);
    tgt = glm::normalize(ders[1]);
}

/**
Evaluate point on a nonrational NURBS surface
@param[in] degree_u Degree of the given surface in u-direction.
@param[in] degree_v Degree of the given surface in v-direction.
@param[in] knots_u Knot vector of the surface in u-direction.
@param[in] knots_v Knot vector of the surface in v-direction.
@param[in] control_points Control points of the surface in a 2d array.
@param[in] u Parameter to evaluate the surface at.
@param[in] v Parameter to evaluate the surface at.
@param[in, out] point Resulting point on the surface at (u, v).
*/
template <int dim, typename T>
void surfacePoint(unsigned int degree_u, unsigned int degree_v,
                  const std::vector<T> &knots_u, const std::vector<T> &knots_v,
                  const array2<glm::vec<dim, T>> &control_points,
                  T u, T v, glm::vec<dim, T> &point) {

    // Initialize result to 0s
    for (int i = 0; i < dim; i++) {
        point[i] = static_cast<T>(0.0);
    }

    // Find span and non-zero basis functions
    int span_u = findSpan(degree_u, knots_u, u);
    int span_v = findSpan(degree_v, knots_v, v);
    std::vector<T> Nu, Nv;
    bsplineBasis(degree_u, span_u, knots_u, u, Nu);
    bsplineBasis(degree_v, span_v, knots_v, v, Nv);

    for (int l = 0; l <= degree_v; l++) {
        glm::vec<dim, T> temp(0.0);
        for (int k = 0; k <= degree_u; k++) {
            temp += static_cast<T>(Nu[k]) *
                    control_points(span_u - degree_u + k, span_v - degree_v + l);
        }

        point += static_cast<T>(Nv[l]) * temp;
    }
}

/**
Evaluate point on a nonrational NURBS surface
@param[in] srf Surface object
@param[in] u Parameter to evaluate the surface at.
@param[in] v Parameter to evaluate the surface at.
@param[in, out] point Resulting point on the surface at (u, v).
*/
template <int dim, typename T>
void surfacePoint(const Surface<dim, T> &srf, T u, T v, glm::vec<dim, T> &point) {
    surfacePoint(srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v, srf.control_points,
                 u, v, point);
}

/**
Evaluate point on a non-rational NURBS surface
@param[in] degree_u Degree of the given surface in u-direction.
@param[in] degree_v Degree of the given surface.
@param[in] knots Knot vector of the surface.
@param[in] control_points Control points of the surface in a 2D array.
@param[in] weights Weights corresponding to each control point.
@param[in] u Parameter to evaluate the surface at.
@param[in] v Parameter to evaluate the surface at.
@param[in, out] point Resulting point on the surface at (u, v).
*/
template <int dim, typename T>
void rationalSurfacePoint(unsigned int degree_u, unsigned int degree_v,
                          const std::vector<T> &knots_u, const std::vector<T> &knots_v,
                          const array2<glm::vec<dim, T>> &control_points,
                          const array2<T> &weights, T u, T v,
                          glm::vec<dim, T> &point) {

    typedef glm::vec<dim + 1, T> tvecnp1;

    // Compute homogenous coordinates of control points
    array2<tvecnp1> Cw;
    Cw.resize(control_points.rows(), control_points.cols());
    for (int i = 0; i < control_points.rows(); i++) {
        for (int j = 0; j < control_points.cols(); j++) {
            Cw(i, j) = tvecnp1(util::cartesianToHomogenous(control_points(i, j),
                                                           weights(i, j)));
        }
    }

    // Compute point using homogenous coordinates
    tvecnp1 pointw;
    surfacePoint(degree_u, degree_v, knots_u, knots_v, Cw, u, v, pointw);

    // Convert back to cartesian coordinates
    point = util::homogenousToCartesian(pointw);
}

/**
Evaluate point on a non-rational NURBS surface
@param[in] srf RationalSurface object
@param[in] u Parameter to evaluate the surface at.
@param[in] v Parameter to evaluate the surface at.
@param[in, out] point Resulting point on the surface at (u, v).
*/
template <int dim, typename T>
void rationalSurfacePoint(const RationalSurface<dim, T> &srf, T u, T v, glm::vec<dim, T> &point) {
    rationalSurfacePoint(srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v, srf.control_points,
                         srf.weights, u, v, point);
}

/**
Evaluate derivatives on a non-rational NURBS surface
@param[in] degree_u Degree of the given surface in u-direction.
@param[in] degree_v Degree of the given surface in v-direction.
@param[in] knots_u Knot vector of the surface in u-direction.
@param[in] knots_v Knot vector of the surface in v-direction.
@param[in] control_points Control points of the surface in a 2D array.
@param[in] num_ders Number of times to differentiate
@param[in] u Parameter to evaluate the surface at.
@param[in] v Parameter to evaluate the surface at.
@param[in, out] surf_ders Derivatives of the surface at (u, v).
*/
template <int dim, typename T>
void surfaceDerivatives(unsigned int degree_u, unsigned int degree_v,
                        const std::vector<T> &knots_u, const std::vector<T> &knots_v,
                        const array2<glm::vec<dim, T>> &control_points, int num_ders,
                        T u, T v, array2<glm::vec<dim, T>> &surf_ders) {

    surf_ders.clear();
    surf_ders.resize(num_ders + 1, num_ders + 1);

    // Set higher order derivatives to 0
    for (int k = degree_u + 1; k <= num_ders; k++) {
        for (int l = degree_v + 1; l <= num_ders; l++) {
            surf_ders(k, l) = glm::vec<dim, T>(0.0);
        }
    }

    // Find span and basis function derivatives
    int span_u = findSpan(degree_u, knots_u, u);
    int span_v = findSpan(degree_v, knots_v, v);
    array2<T> ders_u, ders_v;
    bsplineDerBasis(degree_u, span_u, knots_u, u, num_ders, ders_u);
    bsplineDerBasis(degree_v, span_v, knots_v, v, num_ders, ders_v);

    // Number of non-zero derivatives is <= degree
    int du = num_ders < degree_u ? num_ders : degree_u;
    int dv = num_ders < degree_v ? num_ders : degree_v;
    std::vector<glm::vec<dim, T>> temp;
    temp.resize(degree_v + 1);

    // Compute derivatives
    for (int k = 0; k <= du; k++) {
        for (int s = 0; s <= degree_v; s++) {
            temp[s] = glm::vec<dim, T>(0.0);
            for (int r = 0; r <= degree_u; r++) {
                temp[s] += static_cast<T>(ders_u(k, r)) *
                           control_points(span_u - degree_u + r, span_v - degree_v + s);
            }
        }

        int nk = num_ders - k;
        int dd = nk < dv ? nk : dv;

        for (int l = 0; l <= dd; l++) {
            surf_ders(k, l) = glm::vec<dim, T>(0.0);

            for (int s = 0; s <= degree_v; s++) {
                surf_ders(k, l) += static_cast<T>(ders_v(l, s)) * temp[s];
            }
        }
    }
}

/**
Evaluate derivatives on a non-rational NURBS surface
@param[in] degree_u Degree of the given surface in u-direction.
@param[in] degree_v Degree of the given surface in v-direction.
@param[in] knots_u Knot vector of the surface in u-direction.
@param[in] knots_v Knot vector of the surface in v-direction.
@param[in] control_points Control points of the surface in a 2D array.
@param[in] num_ders Number of times to differentiate
@param[in] u Parameter to evaluate the surface at.
@param[in] v Parameter to evaluate the surface at.
@param[in, out] surf_ders Derivatives of the surface at (u, v).
*/
template <int dim, typename T>
void surfaceDerivatives(const Surface<dim, T> &srf, int num_ders, T u, T v, array2<glm::vec<dim, T>> &surf_ders) {
    surfaceDerivatives(srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v, srf.control_points,
                       srf.weights, num_ders, u, v, surf_ders);
}

/**
Evaluate derivatives on a rational NURBS surface
@param[in] degree_u Degree of the given surface in u-direction.
@param[in] degree_v Degree of the given surface in v-direction.
@param[in] knots_u Knot vector of the surface in u-direction.
@param[in] knots_v Knot vector of the surface in v-direction.
@param[in] control_points Control points of the surface.
@param[in] weights Weights corresponding to each control point.
@param[in] u Parameter to evaluate the surface at.
@param[in] v Parameter to evaluate the surface at.
@param[in] num_ders Number of times to differentiate
@param[in, out] surf_ders Derivatives on the surface at parameter (u, v).
*/
template <int dim, typename T>
void rationalSurfaceDerivatives(unsigned int degree_u, unsigned int degree_v,
                                const std::vector<T> &knots_u, const std::vector<T> &knots_v,
                                const array2<glm::vec<dim, T>> &control_points,
                                const array2<T> &weights, int num_ders, T u, T v,
                                array2<glm::vec<dim, T>> &surf_ders) {

    using namespace std;
    using namespace glm;

    typedef vec<dim, T> tvecn;
    typedef vec<dim + 1, T> tvecnp1;

    array2<tvecnp1> homo_cp;
    homo_cp.resize(control_points.rows(), control_points.cols());
    for (int i = 0; i < control_points.rows(); ++i) {
        for (int j = 0; j < control_points.cols(); ++j) {
            homo_cp(i, j) = util::cartesianToHomogenous(control_points(i, j), weights(i, j));
        }
    }

    array2<tvecnp1> homo_ders;
    surfaceDerivatives(u, v, degree_u, degree_v, knots_u, knots_v, homo_cp, num_ders, homo_ders);

    array2<tvecn> Aders;
    Aders.resize(num_ders + 1, num_ders + 1);
    for (int i = 0; i < homo_ders.rows(); ++i) {
        for (int j = 0; j < homo_ders.cols(); ++j) {
            Aders(i, j) = util::truncateHomogenous(homo_ders(i, j));
        }
    }

    surf_ders.resize(num_ders + 1, num_ders + 1);
    for (int k = 0; k < num_ders + 1; ++k) {
        for (int l = 0; l < num_ders - k + 1; ++l) {
            auto der = Aders(k, l);

            for (int j = 1; j < l + 1; ++j) {
                der -= (T)util::binomial(l, j) * homo_ders(0, j)[dim] * surf_ders(k, l - j);
            }

            for (int i = 1; i <  k + 1; ++i) {
                der -= (T)util::binomial(k, i) * homo_ders(i, 0)[dim] * surf_ders(k - i, l);

                tvecn tmp((T)0.0);
                for (int j = 1; j < l + 1; ++j) {
                    tmp -= (T)util::binomial(l, j) * homo_ders(i, j)[dim] * surf_ders(k - 1, l - j);
                }

                der -= (T)util::binomial(k, i) * tmp;
            }

            der *= 1 / homo_ders(0, 0)[dim];
            surf_ders(k, l) = der;
        }
    }
}

/**
Evaluate derivatives on a rational NURBS surface
@param[in] srf RationalSurface object
@param[in] u Parameter to evaluate the surface at.
@param[in] v Parameter to evaluate the surface at.
@param[in] num_ders Number of times to differentiate
@param[in, out] surf_ders Derivatives on the surface at parameter (u, v).
*/
template <int dim, typename T>
void rationalSurfaceDerivatives(const RationalSurface<dim, T> &srf, int num_ders, T u, T v,
                                array2<glm::vec<dim, T>> &surf_ders) {
    rationalSurfaceDerivatives(srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v, srf.control_points,
                               srf.weights, num_ders, u, v, surf_ders);
}

/**
Evaluate the two orthogonal tangents of a non-rational surface at the given parameters
@param[in] srf Surface object
@param u Parameter in the u-direction
@param v Parameter in the v-direction
@param[in, out] du Unit tangent along u-direction
@param[in, out] dv Unit tangent along v-direction
*/
template <int dim, typename T>
void surfaceTangent(const Surface<dim, T> &srf, T u, T v,
                    glm::vec<dim, T> &du, glm::vec<dim, T> &dv) {
    array2<glm::vec<dim, T>> ptder;
    surfaceDerivatives(srf, 1, u, v, ptder);
    du = ptder[1][0];
    dv = ptder[0][1];
    du = glm::normalize(du);
    dv = glm::normalize(dv);
}

/**
Evaluate the two orthogonal tangents of a rational surface at the given parameters
@param[in] srf Rational Surface object
@param u Parameter in the u-direction
@param v Parameter in the v-direction
@param[in, out] du Unit tangent along u-direction
@param[in, out] dv Unit tangent along v-direction
*/
template <int dim, typename T>
void rationalSurfaceTangent(const RationalSurface<dim, T> &srf, T u, T v,
                            glm::vec<dim, T> &du, glm::vec<dim, T> &dv) {
    array2<glm::vec<dim, T>> ptder;
    rationalSurfaceDerivatives(srf, 1, u, v, ptder);
    du = ptder[1][0];
    dv = ptder[0][1];
    du = glm::normalize(du);
    dv = glm::normalize(dv);
}


/**
Evaluate the normal a non-rational surface at the given parameters
@param[in] srf Surface object
@param u Parameter in the u-direction
@param v Parameter in the v-direction
@param[in, out] normal Unit normal at of the surface at (u, v)
*/
template <int dim, typename T>
void surfaceNormal(const Surface<dim, T> &srf, T u, T v, glm::vec<dim, T> &normal) {
    array2<glm::vec<dim, T>> ptder;
    surfaceDerivatives(srf, 1, u, v, ptder);
    normal = glm::cross(ptder[0][1], ptder[1][0]);
    normal = glm::normalize(normal);
}

/**
Evaluate the normal of a rational surface at the given parameters
@param[in] srf Rational Surface object
@param u Parameter in the u-direction
@param v Parameter in the v-direction
@param[in, out] normal Unit normal at of the surface at (u, v)
*/
template <int dim, typename T>
void rationalSurfaceNormal(const RationalSurface<dim, T> &srf, T u, T v, glm::vec<dim, T> &normal) {
    array2<glm::vec<dim, T>> ptder;
    rationalSurfaceDerivatives(srf, 1, u, v, ptder);
    normal = glm::cross(ptder[0][1], ptder[1][0]);
    normal = glm::normalize(normal);
}

} // namespace nurbstk
