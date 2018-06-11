/*
@file nurbstk/surface.h
@author Pradeep Kumar Jayaraman <pradeep.pyro@gmail.com>

The Surface and RationalSurface classes represent non-rational and rational NURBS surfaces, respectively.

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE.txt file.
*/

#pragma once
#include <vector>
#include <stdexcept>

#include "evaluate.h"
#include "knots.h"
#include "io.h"

namespace nurbstk {

/**
\brief Base Class for representing a NURBS surface
\tparam nd Dimension of the surface (always 3)
\tparam T Data type of control points and weights (float or double)
*/
template <int nd, typename T>
class BaseSurface {
public:
    typedef glm::vec<nd, T> vecnt;

    explicit BaseSurface(unsigned int degreeU, unsigned int degreeV,
                         const std::vector<T> &knotsU, const std::vector<T> &knotsV,
                         const std::vector<std::vector<vecnt>> &controlPoints) {

        if (degreeU < 1 && degreeV < 1) {
            throw std::logic_error("Degree has to be atleast 1");
        }

        if (!isValidRelation(degreeU, knotsU.size(), controlPoints.size()) ||
                !isValidRelation(degreeV, knotsV.size(), controlPoints[0].size())) {
            throw std::logic_error("Invalid relation between degree, "
                                   "knots and control points");
        }

        if (!isKnotVectorMonotonic(knotsU) || !isKnotVectorMonotonic(knotsV)) {
            throw(std::logic_error("Knot vector(s) is not monotonic"));
        }

        this->degree_u_ = degreeU;
        this->degree_v_ = degreeV;
        this->knots_u_ = knotsU;
        this->knots_v_ = knotsV;
        this->control_points_ = controlPoints;
    }

    /**
    Returns the point on the surface at the given parameters
    @param u Parameter in the u-direction
    @param v Parameter in the v-direction
    @return Point on surface at (u, v)
    */
    virtual vecnt point(T u, T v) const = 0;

    /**
    Computes all nDers derivatives of the surface at the given parameters
    @param u Parameter in the u-direction
    @param v Parameter in the v-direction
    @param nDers Number of derivatives to compute
    E.g. nDers=2 for upto 2nd derivative
    @return ders Derivatives at (u, v)
    ders[k][l] is the derivative k times w.r.t. u and l times w.r.t. v
    E.g. ders[0][0] is the point on the surface
    ders[0][1] is the first partial derivative w.r.t. v
    ders[2][0] is the second partial derivative w.r.t. u
    */
    virtual void derivatives(T u, T v, int nDers,
                             std::vector<std::vector<vecnt>> &ders) const = 0;

    /**
    Returns the two orthogonal tangents at the given parameters
    @param u Parameter in the u-direction
    @param v Parameter in the v-direction
    @param normalize Whether to return a unit vector
    */
    void tangent(T u, T v, vecnt &du, vecnt &dv,
                 bool normalize = true) const {
        std::vector<std::vector<vecnt>> ptder;
        derivatives(u, v, 1, ptder);
        du = ptder[1][0];
        dv = ptder[0][1];
        if (normalize) {
            glm::normalize(du);
            glm::normalize(dv);
        }
    }

    /**
    Returns the normal at the given parameters
    @param u Parameter in the u-direction
    @param v Parameter in the v-direction
    @param normalize Whether to return a unit vector
    */
    vecnt normal(T u, T v, bool normalize = true) const {
        std::vector<std::vector<vecnt>> ptder;
        derivatives(u, v, 1, ptder);
        vecnt normal = glm::cross(ptder[0][1], ptder[1][0]);
        return normalize ? glm::normalize(normal) : normal;
    }

    /**
    Returns a 2D vector of control points
    */
    const std::vector<std::vector<vecnt>> & controlPoints() const {
        return control_points_;
    }

    /**
    Returns the degree along u-direction
    */
    unsigned int degreeU() const {
        return degree_u_;
    }

    /**
    Returns the degree along v-direction
    */
    unsigned int degreeV() const {
        return degree_v_;
    }

    /**
    Returns a copy of the knot vector along u-direction
    */
    const std::vector<double> & knotsU() const {
        return knots_u_;
    }

    /**
    Returns knot value at index along u-direction
    */
    T knotU(size_t index) const {
        return knots_u_[index];
    }

    /**
    Sets knot value at index along u-direction
    */
    T & knotU(size_t index) {
        return knots_u_[index];
    }

    /**
    Returns knot value at index along v-direction
    */
    T knotV(size_t index) const {
        return knots_v_[index];
    }

    /**
    Sets knot value at index along v-direction
    */
    T & knotV(size_t index, double val) {
        return knots_v_[index];
    }

    /**
    Returns a copy of the knot vector along v-direction
    */
    const std::vector<double> & knotsV() const {
        return knots_v_;
    }

    /**
    Returns number of knots along u-direction
    */
    size_t numKnotsU() const {
        return knots_u_.size();
    }

    /**
    Returns number of knots along u - direction
    */
    size_t numKnotsV() const {
        return knots_v_.size();
    }

    /**
    Returns number of control points along u-direction
    */
    size_t numControlPointsU() const {
        return control_points_.size();
    }

    /**
    Returns number of control points along v-direction
    */
    size_t numControlPointsV() const {
        if (control_points_.empty()) {
            return 0;
        }
        return control_points_[0].size();
    }

    /**
    Returns the control point at (i, j)
    @param i Index along u-direction
    @param j Index along v-direction
    */
    vecnt controlPoint(size_t i, size_t j) const {
        return control_points_[i][j];
    }

    /**
    Sets the control point at (i, j)
    @param i Index along u-direction
    @param j Index along v-direction
    @param pt New control point at (i, j)
    */
    vecnt & controlPoint(size_t i, size_t j) {
        return control_points_[i][j];
    }

    /**
    Returns whether the surface is closed along the u-direction
    */
    bool isClosedU() const {
        for (int j = 0; j < numControlPointsV(); j++) {
            if (control_points_[0][j] != control_points_[numControlPointsU() - 1][j]) {
                return false;
            }
        }

        if (!isKnotVectorClosed(degree_u_, knots_u_)) {
            return false;
        }
        return true;
    }

    /**
    Returns whether the surface is closed along the v-direction
    */
    bool isClosedV() const {
        for (int i = 0; i < numControlPointsU(); i++) {
            if (control_points_[i][0] != control_points_[i][numControlPointsV() - 1]) {
                return false;
            }
        }

        if (!isKnotVectorClosed(degree_v_, knots_v_)) {
            return false;
        }
        return true;
    }

    /**
    Returns whether the surface is valid by checking the relationship
    between the degree, knots and control points as well as the monotonicity
    of the knot vectors
    */
    bool isValid() const {
        return isValidRelation(degree_u_, knots_u_.size(), control_points_.size()) &&
               isValidRelation(degree_v_, knots_v_.size(), control_points_[0].size()) &&
               isKnotVectorMonotonic(knots_u_) && isKnotVectorMonotonic(knots_v_);
    }

protected:
    unsigned int degree_u_, degree_v_;
    std::vector<T> knots_u_, knots_v_;
    std::vector<std::vector<vecnt>> control_points_;
};

/**
\brief Class for representing a non-rational B-spline surface
\tparam nd Dimension of the surface (always 3)
\tparam T Data type of control points and weights (float or double)
*/
template <int nd, typename T>
class Surface : public BaseSurface<nd, T> {
public:
    typedef glm::vec<nd, T> vecnt;

    /**
    Create a non-rational NURBS surface
    @param degreeU Degree of the surface in u-direction
    @param degreeV Degree of the surface in v-direction
    @param knotsU Knot vector in u-direction
    @param knotsV Knot vector in v-direction
    @param controlPoints 2D vector of control points
    */
    Surface(unsigned int degreeU, unsigned int degreeV,
            const std::vector<T> &knotsU, const std::vector<T> &knotsV,
            const std::vector<std::vector<vecnt>> &controlPoints)
        : BaseSurface<nd, T>(degreeU, degreeV, knotsU, knotsV, controlPoints) {
    }

    static Surface<nd, T> fromOBJ(const std::string &filename) {
        unsigned int deg_u, deg_v;
        std::vector<double> knots_u, knots_v;
        std::vector<std::vector<glm::dvec3>> ctrl_pts;
        std::vector<std::vector<double>> weights;
        bool rational;
        readOBJ(filename, deg_u, deg_v, knots_u, knots_v, ctrl_pts, weights, rational);
        return Surface<nd, T>(deg_u, deg_v, knots_u, knots_v, ctrl_pts);
    }

    vecnt point(T u, T v) const override {
        vecnt pt;
        surfacePoint<nd, T>(u, v, this->degree_u_, this->degree_v_, this->knots_u_, this->knots_v_,
                            this->control_points_, pt);
        return pt;
    }

    void derivatives(T u, T v, int nDers,
                     std::vector<std::vector<vecnt>> &ders) const override {
        surfaceDerivatives<nd, T>(u, v, this->degree_u_, this->degree_v_,
                                  this->knots_u_, this->knots_v_, this->control_points_, nDers, ders);
    }
};

/**
\brief Class for representing a rational B-spline (NURBS) surface
\tparam nd Dimension of the surface (always 3)
\tparam T Data type of control points and weights (float or double)
*/
template <int nd, typename T>
class RationalSurface : public BaseSurface<nd, T> {
public:
    typedef glm::vec<nd, T> vecnt;

    /**
    Create a rational NURBS surface
    @param degreeU Degree of the surface in u-direction
    @param degreeV Degree of the surface in v-direction
    @param knotsU Knot vector in u-direction
    @param knotsV Knot vector in v-direction
    @param controlPoints 2D vector of control points
    @param weights 2D vector of weights corresponding to control points
    */
    RationalSurface(unsigned int degreeU, unsigned int degreeV,
                    const std::vector<T> &knotsU, const std::vector<T> &knotsV,
                    const std::vector<std::vector<vecnt>> &controlPoints,
                    const std::vector<std::vector<T>> &weights)
        : BaseSurface<nd, T>(degreeU, degreeV, knotsU, knotsV, controlPoints) {
        if (weights.size() != controlPoints.size()) {
            throw std::logic_error("Number of weights should be equal to number of control points");
        }
        this->weights_ = weights;
    }

    static RationalSurface<nd, T> fromOBJ(const std::string &filename) {
        unsigned int deg_u, deg_v;
        std::vector<double> knots_u, knots_v;
        std::vector<std::vector<glm::dvec3>> ctrl_pts;
        std::vector<std::vector<double>> weights;
        bool rational;
        readOBJ(filename, deg_u, deg_v, knots_u, knots_v, ctrl_pts, weights, rational);
        return RationalSurface<nd, T>(deg_u, deg_v, knots_u, knots_v, ctrl_pts, weights);
    }

    /**
    Returns the weight of the control point at (i, j)
    @param i Index along u-direction
    @param j Index along v-direction
    */
    T weight(size_t i, size_t j) const {
        return weights_[i][j];
    }

    /**
    Returns a mutable reference to the weight of the control point at (i, j)
    @param i Index along u-direction
    @param j Index along v-direction
    */
    T & weight(size_t i, size_t j) {
        return weights_[i][j];
    }

    vecnt point(T u, T v) const override {
        vecnt pt;
        rationalSurfacePoint<nd, T>(u, v, this->degree_u_, this->degree_v_, this->knots_u_, this->knots_v_,
                                    this->control_points_, weights_, pt);
        return pt;
    }

    void derivatives(T u, T v, int nDers, std::vector<std::vector<vecnt>> &ders) const override {
        rationalSurfaceDerivatives<nd, T>(u, v, this->degree_u_, this->degree_v_, this->knots_u_, this->knots_v_,
                                          this->control_points_, weights_, nDers, ders);
    }

private:
    std::vector<std::vector<T>> weights_;
};

// Typedefs for ease of use
typedef Surface<3, float> Surface3f;
typedef Surface<3, double> Surface3d;
typedef RationalSurface<3, float> RationalSurface3f;
typedef RationalSurface<3, double> RationalSurface3d;

} // namespace nurbstk
