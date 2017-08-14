/*
@file openspline/nurbs_surface.h
@author Pradeep Kumar Jayaraman <pradeep.pyro@gmail.com>

The NurbsSurface class represents an non-uniform rational B-spline surface.

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE.txt file.
*/

#pragma once
#include <vector>
#include <exception>

#include "surface.h"
#include "nurbs_evaluate.h"
#include "nurbs_knots.h"

namespace ospl {

/**
\brief Class for representing a NURBS surface
\tparam nd Dimension of the surface
\tparam T Data type of control points and weights (float or double)
*/
template <int nd, typename T>
class NurbsSurface : public Surface<nd, T> {
public:
	typedef glm::vec<nd, T> vecnt;

	/**
	Create a non-rational NurbsSurface
	@param degreeU Degree of the surface in u-direction
	@param degreeV Degree of the surface in v-direction
	@param knotsU Knot vector in u-direction
	@param knotsV Knot vector in v-direction
	@param controlPoints 2D vector of control points
	*/
	explicit NurbsSurface(unsigned int degreeU, unsigned int degreeV,
		const std::vector<double> &knotsU, const std::vector<double> &knotsV,
		const std::vector<std::vector<vecnt>> &controlPoints) {

		if (degreeU < 1 && degreeV < 1) {
			throw std::logic_error("Degree has to be atleast 1");
		}

		if (!isValidRelation(degreeU, knotsU.size(), controlPoints.size()) ||
			!isValidRelation(degreeV, knotsV.size(), controlPoints[0].size())) {
			throw std::logic_error("Invalid relation between degree, "
				"knots and control points");
		}

		this->degU = degreeU;
		this->degV = degreeV;
		this->knotsU = knotsU;
		this->knotsV = knotsV;
		this->cp = controlPoints;
		this->isRat = false;

		if (this->w.empty()) {
			this->w.resize(this->cp.size());
			for (auto &vec : this->w) {
				vec.resize(this->cp[0].size(), 1.0);
			}
		}
	}

	/**
	Create a rational NurbsSurface
	@param degreeU Degree of the surface in u-direction
	@param degreeV Degree of the surface in v-direction
	@param knotsU Knot vector in u-direction
	@param knotsV Knot vector in v-direction
	@param controlPoints 2D vector of control points
	@param weights 2D vector of weights corresponding to control points
	*/
	explicit NurbsSurface(unsigned int degreeU, unsigned int degreeV,
		const std::vector<double> &knotsU, const std::vector<double> &knotsV,
		const std::vector<std::vector<vecnt>> &controlPoints,
		const std::vector<std::vector<T>> &weights) {

		if (degreeU < 1 && degreeV < 1) {
			throw std::logic_error("Degree has to be atleast 1");
		}

		if (!isValidRelation(degreeU, knotsU.size(), controlPoints.size()) ||
			!isValidRelation(degreeV, knotsV.size(), controlPoints[0].size())) {
			throw std::logic_error("Invalid relation between degree, "
				"knots and control points");
		}

		this->degU = degreeU;
		this->degV = degreeV;
		this->knotsU = knotsU;
		this->knotsV = knotsV;
		this->cp = controlPoints;
		this->isRat = true;

		this->w = weights;
	}

	/**
	Returns the point on the surface at the given parameters
	@param u Parameter in the u-direction
	@param v Parameter in the v-direction
	@return Point on surface at (u, v)
	*/
	vecnt point(double u, double v) const override {
		vecnt pt;
		if (!rational()) {
			nurbsSurfacePoint<nd, T>(u, v, degU, degV, knotsU, knotsV, cp, pt);
		}
		else {
			nurbsRationalSurfacePoint<nd, T>(u, v, degU, degV, knotsU, knotsV,
				cp, w, pt);
		}
		return pt;
	}

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
	void derivatives(double u, double v, int nDers,
		std::vector<std::vector<vecnt>> &ders) const override {
		if (!rational()) {
			nurbsSurfaceDerivatives<nd, T>(u, v, degU, degV,
				knotsU, knotsV, cp, nDers, ders);
		}
		else {
			nurbsRationalSurfaceDerivatives<nd, T>(u, v, degU, degV, knotsU,
				knotsV, cp, w, nDers, ders);
		}
	}

	/**
	Returns the two orthogonal tangents at the given parameters
	@param u Parameter in the u-direction
	@param v Parameter in the v-direction
	@param normalize Whether to return a unit vector
	*/
	void tangent(double u, double v, vecnt &du, vecnt &dv,
		bool normalize = true) const override {
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
	vecnt normal(double u, double v, bool normalize = true) const override {
		std::vector<std::vector<vecnt>> ptder;
		derivatives(u, v, 1, ptder);
		vecnt normal = glm::cross(ptder[0][1], ptder[1][0]);
		return normalize ? glm::normalize(normal) : normal;
	}

	/**
	Returns a 2D vector of control points
	*/
	const std::vector<std::vector<vecnt>> & controlPoints() const {
		return cp;
	}

	/**
	Returns the degree along u-direction
	*/
	unsigned int degreeU() const override {
		return degU;
	}

	/**
	Returns the degree along v-direction
	*/
	unsigned int degreeV() const override {
		return degV;
	}

	/**
	Returns knot vector along u-direction
	*/
	const std::vector<double> & knotVectorU() const {
		return knotsU;
	}

	/**
	Returns knot vector along v-direction
	*/
	const std::vector<double> & knotVectorV() const {
		return knotsV;
	}

	/**
	Returns number of control points along u-direction
	*/
	size_t numControlPointsU() const {
		return cp.size();
	}

	/**
	Returns number of control points along v-direction
	*/
	size_t numControlPointsV() const {
		if (cp.empty()) {
			return 0;
		}
		return cp[0].size();
	}

	/**
	Returns the control point at (i, j)
	@param i Index along u-direction
	@param j Index along v-direction
	*/
	vecnt controlPoint(size_t i, size_t j) const {
		return cp[i][j];
	}
	
	/**
	Sets the control point at (i, j)
	@param i Index along u-direction
	@param j Index along v-direction
	@param pt New control point at (i, j)
	*/
	void setControlPoint(size_t i, size_t j, vecnt pt) {
		cp[i][j] = pt;
	}

	T weight(size_t i, size_t j) const {
		return w[i][j];
	}

	void setWeight(size_t i, size_t j, T weight) {
		w[i][j] = weight;
	}

	/**
	Returns whether the surface is rational
	*/
	bool rational() const {
		return isRat;
	}

	/**
	Make the surface rational/non-rational
	*/
	void setRational(bool flag) {
		isRat = flag;
	}

	/**
	Returns whether the surface is closed along the u-direction
	*/
	bool isClosedU() const {
		for (int j = 0; j < numControlPointsV(); j++) {
			if (cp[0][j] != cp[numControlPointsU() - 1][j]) {
				return false;
			}
		}
		
		if (!isKnotVectorClosed(degU, knotsU)) {
			return false;
		}
		return true;
	}

	/**
	Returns whether the surface is closed along the v-direction
	*/
	bool isClosedV() const {
		for (int i = 0; i < numControlPointsU(); i++) {
			if (cp[i][0] != cp[i][numControlPointsV() - 1]) {
				return false;
			}
		}

		if (!isKnotVectorClosed(degV, knotsV)) {
			return false;
		}
		return true;
	}

private:
	unsigned int degU, degV;
	std::vector<double> knotsU, knotsV;
	std::vector<std::vector<vecnt>> cp;
	std::vector<std::vector<T>> w;
	bool isRat;
};

// Typedefs for ease of use
typedef NurbsSurface<3, float> NurbsSurface3f;
typedef NurbsSurface<3, double> NurbsSurface3d;

} // namespace ospl