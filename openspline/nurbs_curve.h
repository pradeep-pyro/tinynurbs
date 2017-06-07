/*
@file openspline/nurbs_curve.h
@author Pradeep Kumar Jayaraman <pradeep.pyro@gmail.com>

The NurbsCurve class represents an non-uniform rational B-spline curve.

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE.txt file.
*/

#pragma once
#include <vector>
#include <exception>

#include "nurbs_evaluate.h"

namespace ospl {

template <int nd, typename T>
class NurbsCurve {
public:
	typedef glm::vec<nd, T> vecnt;

	NurbsCurve(unsigned int degree, const std::vector<double> &knots,
		const std::vector<vecnt> &controlPoints) {
		if (!isValidRelation(degree, knots.size(), 
			controlPoints.size())) {
			throw std::logic_error("nKnots != degree + nCtrlPts + 1");
		}
		this->deg = degree;
		this->knots = knots;
		this->cp = controlPoints;
		this->w.resize(controlPoints.size(), 1.0);
		this->isRat = false;
	}

	NurbsCurve(unsigned int degree, const std::vector<double> &knots, 
		const std::vector<vecnt> &controlPoints, const std::vector<T> &weights) {
		NurbsCurve(degree, knots, controlPoints);
		if (!isValidRelation(degree, knots.size(),
			controlPoints.size())) {
			throw std::logic_error("nKnots != degree + nCtrlPts + 1");
		}
		this->deg = degree;
		this->knots = knots;
		this->cp = controlPoints;
		this->w = weights;
		this->isRat = true;
	}

	int degree() const {
		return deg;
	}

	vecnt point(double u) {
		vecnt pt;
		if (!rational()) {
			nurbsCurvePoint<nd, T>(u, deg, knots, cp, pt);
		}
		else {
			nurbsRationalCurvePoint<nd, T>(u, deg, knots, cp, w, pt);
		}
		return pt;
	}

	void derivatives(double u, int nDers, std::vector<vecnt> &ptder) {
		if (!rational()) {
			nurbsCurveDerivatives<nd, T>(u, deg, knots, cp, nDers, ptder);
		}
		else {
			nurbsRationalCurveDerivatives<nd, T>(u, deg, knots, cp, w, 1, ptder);
		}
	}

	vecnt tangent(double u, bool normalize = true) {
		std::vector<vecnt> ders;
		derivatives(u, 1, ders);
		return normalize ? glm::normalize(ders[1]) : ders[1];
	}

	vecnt controlPoint(int i) const {
		return cp[i];
	}

	void setControlPoint(int i, vecnt pt) {
		cp[i] = pt;
	}

	const std::vector<vecnt>& controlPoints() const {
		return cp;
	}

	size_t numControlPoints() const {
		return cp.size();
	}

	T weight(int i) const {
		return w[i];
	}

	void setWeight(int i, T weight) {
		w[i] = weight;
	}

	const std::vector<double>& knotVector() const {
		return knots;
	}

	bool rational() const {
		return isRat;
	}

	void setRational(bool isRational) {
		isRat = isRational;
	}
	
private:
	unsigned int deg;
	std::vector<double> knots;
	std::vector<vecnt> cp;
	std::vector<T> w;
	bool isRat;
};

// Typedefs for ease of use
typedef NurbsCurve<2, float> NurbsCurve2f;
typedef NurbsCurve<2, double> NurbsCurve2d;
typedef NurbsCurve<3, float> NurbsCurve3f;
typedef NurbsCurve<3, double> NurbsCurve3d;

} // namespace ospl