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

namespace nurbstk {

/**
\brief Base Class for representing a NURBS surface
\tparam nd Dimension of the surface
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

        this->deg_u = degreeU;
        this->deg_v = degreeV;
        this->knots_u = knotsU;
        this->knots_v = knotsV;
        this->ctrlPts = controlPoints;
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
        return ctrlPts;
    }

    /**
    Returns the degree along u-direction
    */
    unsigned int degreeU() const {
        return deg_u;
    }

    /**
    Returns the degree along v-direction
    */
    unsigned int degreeV() const {
        return deg_v;
    }

    /**
    Returns a copy of the knot vector along u-direction
    */
    const std::vector<double> & knotsU() const {
        return knots_u;
    }

    /**
    Returns knot value at index along u-direction
    */
    T knotU(size_t index) const {
        return knots_u[index];
    }

    /**
    Sets knot value at index along u-direction
    */
    T & knotU(size_t index) {
        return knots_u[index];
    }

    /**
    Returns knot value at index along v-direction
    */
    T knotV(size_t index) const {
        return knots_v[index];
    }

    /**
    Sets knot value at index along v-direction
    */
    T & knotV(size_t index, double val) {
        return knots_v[index];
    }

    /**
    Returns a copy of the knot vector along v-direction
    */
    const std::vector<double> & knotsV() const {
        return knots_v;
    }

    /**
    Returns number of knots along u-direction
    */
    size_t numKnotsU() const {
        return knots_u.size();
    }

    /**
    Returns number of knots along u - direction
    */
    size_t numKnotsV() const {
        return knots_v.size();
    }

    /**
    Returns number of control points along u-direction
    */
    size_t numControlPointsU() const {
        return ctrlPts.size();
    }

    /**
    Returns number of control points along v-direction
    */
    size_t numControlPointsV() const {
        if (ctrlPts.empty()) {
            return 0;
        }
        return ctrlPts[0].size();
    }

    /**
    Returns the control point at (i, j)
    @param i Index along u-direction
    @param j Index along v-direction
    */
    vecnt controlPoint(size_t i, size_t j) const {
        return ctrlPts[i][j];
    }

    /**
    Sets the control point at (i, j)
    @param i Index along u-direction
    @param j Index along v-direction
    @param pt New control point at (i, j)
    */
    vecnt & controlPoint(size_t i, size_t j) {
        return ctrlPts[i][j];
    }

    /**
    Returns whether the surface is closed along the u-direction
    */
    bool isClosedU() const {
        for (int j = 0; j < numControlPointsV(); j++) {
            if (ctrlPts[0][j] != ctrlPts[numControlPointsU() - 1][j]) {
                return false;
            }
        }

        if (!isKnotVectorClosed(deg_u, knots_u)) {
            return false;
        }
        return true;
    }

    /**
    Returns whether the surface is closed along the v-direction
    */
    bool isClosedV() const {
        for (int i = 0; i < numControlPointsU(); i++) {
            if (ctrlPts[i][0] != ctrlPts[i][numControlPointsV() - 1]) {
                return false;
            }
        }

        if (!isKnotVectorClosed(deg_v, knots_v)) {
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
        return isValidRelation(deg_u, knots_u.size(), ctrlPts.size()) &&
               isValidRelation(deg_v, knots_v.size(), ctrlPts[0].size()) &&
               isKnotVectorMonotonic(knots_u) && isKnotVectorMonotonic(knots_v);
    }

protected:
    unsigned int deg_u, deg_v;
    std::vector<T> knots_u, knots_v;
    std::vector<std::vector<vecnt>> ctrlPts;
};


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

    vecnt point(T u, T v) const override {
        vecnt pt;
        surfacePoint<nd, T>(u, v, this->deg_u, this->deg_v, this->knots_u, this->knots_v,
                            this->ctrlPts, pt);
        return pt;
    }

    void derivatives(T u, T v, int nDers,
                     std::vector<std::vector<vecnt>> &ders) const override {
        surfaceDerivatives<nd, T>(u, v, this->deg_u, this->deg_v,
                                  this->knots_u, this->knots_v, this->ctrlPts, nDers, ders);
    }
};


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
        this->weights = weights;
    }

    /**
    Returns the weight of the control point at (i, j)
    @param i Index along u-direction
    @param j Index along v-direction
    */
    T weight(size_t i, size_t j) const {
        return weights[i][j];
    }

    /**
    Returns a mutable reference to the weight of the control point at (i, j)
    @param i Index along u-direction
    @param j Index along v-direction
    */
    T & weight(size_t i, size_t j) {
        return weights[i][j];
    }

    vecnt point(T u, T v) const override {
        vecnt pt;
        rationalSurfacePoint<nd, T>(u, v, this->deg_u, this->deg_v, this->knots_u, this->knots_v,
                                    this->ctrlPts, weights, pt);
        return pt;
    }

    void derivatives(T u, T v, int nDers, std::vector<std::vector<vecnt>> &ders) const override {
        rationalSurfaceDerivatives<nd, T>(u, v, this->deg_u, this->deg_v, this->knots_u, this->knots_v,
                                          this->ctrlPts, weights, nDers, ders);
    }
private:
    std::vector<std::vector<T>> weights;
};

// Typedefs for ease of use
typedef Surface<3, float> Surface3f;
typedef Surface<3, double> Surface3d;
typedef RationalSurface<3, float> RationalSurface3f;
typedef RationalSurface<3, double> RationalSurface3d;

} // namespace nurbstk
