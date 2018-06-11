/*
@file nurbstk/curve.h
@author Pradeep Kumar Jayaraman <pradeep.pyro@gmail.com>

The Curve class represents a non-uniform polynomial B-spline curve, while the RationalCurve class
represents a non-uniform rational B-spline (NURBS) curve.

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE.txt file.
*/

#pragma once
#include <vector>
#include <exception>
#include <stdexcept>

#include "evaluate.h"
#include "knots.h"

namespace nurbstk {

/**
@brief Base class for representing a polynomial and rational
@tparam nd Dimension of the curve (2 or 3)
@tparam T Data type of control points and knots (float or double)
*/
template <int nd, typename T>
class BaseCurve {
public:
    typedef glm::vec<nd, T> vecnt;

    BaseCurve(unsigned int degree, const std::vector<T> &knots,
              const std::vector<vecnt> &control_points) {
        if (degree < 1) {
            throw std::logic_error("Degree has to be atleast 1");
        }
        if (!isValidRelation(degree, knots.size(), control_points.size())) {
            throw std::logic_error("Invalid relation between degree, "
                                   "knots and control points");
        }
        if (!isKnotVectorMonotonic(knots)) {
            throw std::logic_error("Knot vector has to be monotonic");
        }
        degree_ = degree;
        knots_ = knots;
        control_points_ = control_points;
    }

    /**
    Get the degree of the curve
    */
    unsigned int degree() const {
        return degree_;
    }

    /**
    Get a control point defining the curve
    @param i Index of the control point
    */
    const vecnt & controlPoint(int i) const {
        return control_points_[i];
    }

    /**
    Get a mutable reference to the control point defining the curve
    @param i Index of the control point
    */
    vecnt & controlPoint(int i) {
        return control_points_[i];
    }

    /**
    Get all the control points defining the curve
    @param i Index of the control point
    */
    const std::vector<vecnt>& controlPoints() const {
        return control_points_;
    }

    /**
    Get the number of control points
    */
    size_t numControlPoints() const {
        return control_points_.size();
    }

    /**
    Get a knot value
    @param i Index of the knot in the knot vector
    */
    T knot(int i) const {
        return knots_[i];
    }

    /**
    Get a mutable reference to a knot value
    @param i Index of the knot in the knot vector
    */
    T & knot(int i) {
        return knots_[i];
    }

    /**
    Get all the knots defining the curve
    */
    const std::vector<T>& knots() const {
        return knots_;
    }

    /**
    Get the number of knots
    */
    size_t numKnots() const {
        return knots_.size();
    }

    virtual vecnt point(T u) const = 0;

    virtual void derivatives(T u, int nDers, std::vector<vecnt> &ptder) const = 0;

    /**
    Compute the tangent of the curve
    @param u Parameter to compute the tangent
    @param normalize Whether to normalize the tangent vector
    */
    vecnt tangent(T u, bool normalize) const {
        std::vector<vecnt> ders;
        derivatives(u, 1, ders);
        return normalize ? glm::normalize(ders[1]) : ders[1];
    }

protected:
    unsigned int degree_;
    std::vector<T> knots_;
    std::vector<vecnt> control_points_;
};

/**
@brief Class for representing a non-rational B-spline curve
@tparam nd Dimension of the curve (2 or 3)
@tparam T Data type of control points and knots (float or double)
*/
template <int nd, typename T>
class Curve : public BaseCurve<nd, T> {
public:
    typedef glm::vec<nd, T> vecnt;
    /**
    Create a non-rational B-spline curve
    @param degree Degree of the curve (must be >1)
    @param knots Knot vector (must be monotonic)
    @param controlPoints 2D/3D vector of control points
    */
    Curve(unsigned int degree, const std::vector<T> &knots, const std::vector<vecnt> &control_points)
        : BaseCurve<nd, T>(degree, knots, control_points) {
    }

    /**
    Compute point on the curve
    @param u Parameter to compute the point
    */
    vecnt point(T u) const override {
        vecnt pt;
        curvePoint<nd, T>(u, this->degree_, this->knots_, this->control_points_, pt);
        return pt;
    }

    /**
    Compute the derivatives of the curve
    @param u Parameter to compute the derivatives
    @param num_derivs Number of derivatives to compute
    @param[out] ptder Derivatives at u (the vector index denotes the order of the derivative)
    */
    void derivatives(T u, int nDers, std::vector<vecnt> &ptder) const override {
        curveDerivatives<nd, T>(u, this->degree_, this->knots_, this->control_points_, nDers, ptder);
    }
};

/**
@brief Class for representing a rational B-spline (NURBS) curve
@tparam nd Dimension of the curve (2 or 3)
@tparam T Data type of control points and knots (float or double)
*/
template <int nd, typename T>
class RationalCurve : public BaseCurve<nd, T> {
public:
    typedef glm::vec<nd, T> vecnt;
    /**
    Create a rational B-spline curve
    @param degree Degree of the curve (must be >1)
    @param knots Knot vector (must be monotonic)
    @param controlPoints 2D/3D vector of control points
    */
    RationalCurve(unsigned int degree, const std::vector<T> &knots, const std::vector<vecnt> &control_points,
                  const std::vector<T> &weights)
        : BaseCurve<nd, T>(degree, knots, control_points) {
        if (control_points.size() != weights.size()) {
            throw std::logic_error("Number of weights should be equal to number of control points");
        }
        weights_ = weights;
    }

    /**
    Convert a non-rational B-spline curve into a rational curve with unit weights
    @param curve Non-rational curve
    */
    RationalCurve(const Curve<nd, T> &curve)
        : BaseCurve<nd, T>(curve.degree(), curve.knots(), curve.controlPoints()) {
        std::vector<T> weights(curve.numControlPoints(), 1);
        weights_ = weights;
    }

    /**
    Convert a non-rational B-spline curve into a rational curve with given weights
    @param curve Non-rational curve
    @param weights Weights corresponding to each of the control points of curve
    */
    RationalCurve(const Curve<nd, T> & curve, const std::vector<T> &weights)
        : BaseCurve<nd, T>(curve.degree(), curve.knots(), curve.controlPoints()) {
        if (curve.controlPoints().size() != weights.size()) {
            throw std::logic_error("Number of weights should be equal to number of control points");
        }
        weights_ = weights;
    }

    /**
    Compute point on the rational curve
    @param u Parameter to compute the point
    */
    vecnt point(T u) const override {
        vecnt pt;
        rationalCurvePoint<nd, T>(u, this->degree_, this->knots_, this->control_points_, this->weights_, pt);
        return pt;
    }

    /**
    Compute the derivatives of the rational curve
    @param u Parameter to compute the derivatives
    @param num_derivs Number of derivatives to compute
    @param[out] ptder Derivatives at u (the vector index denotes the order of the derivative)
    */
    void derivatives(T u, int num_deriv, std::vector<vecnt> &ptder) const override {
        rationalCurveDerivatives<nd, T>(u, this->degree_, this->knots_, this->control_points_, this->weights_, num_deriv, ptder);
    }

    /**
    Get weight associated with the control point
    @param i Index of the control point
    */
    T weight(int i) const {
        return weights_[i];
    }

    /**
    Get a mutable reference to the weight associated with the control point
    @param i Index of the control point
    */
    T & weight(int i) {
        return weights_[i];
    }
private:
    std::vector<T> weights_;
};

// Typedefs for ease of use
typedef Curve<2, float> Curve2f;
typedef Curve<2, double> Curve2d;
typedef Curve<3, float> Curve3f;
typedef Curve<3, double> Curve3d;
typedef RationalCurve<2, float> RationalCurve2f;
typedef RationalCurve<2, double> RationalCurve2d;
typedef RationalCurve<3, float> RationalCurve3f;
typedef RationalCurve<3, double> RationalCurve3d;

} // namespace nurbstk
