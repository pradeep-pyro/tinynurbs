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

	NurbsCurve() {};

	NurbsCurve(unsigned int degree, const std::vector<double> &knots,
		const std::vector<vecnt> &controlPoints) {
		if (!nurbs::isValidRelation(degree, knots.size(), 
			controlPoints.size())) {
			throw std::logic_error("nKnots != degree + nCtrlPts + 1");
		}
		deg_ = degree;
		knots_ = knots;
		cp_ = controlPoints;
		isRat_ = false;
	}

	NurbsCurve(unsigned int degree, const std::vector<double> &knots, 
		const std::vector<vecnt> &controlPoints, const std::vector<T> &weights)
		throw (std::logic_error) {
		NurbsCurve(degree, knots, controlPoints);
		w_ = weights;
		isRat_ = true;
	}
	
	vecnt point(double u) {
		vecnt pt;
		if (!isRat_) {
			nurbs::curvePoint<nd, T>(u, deg_, knots_, cp_, pt);
		}
		else {
			nurbs::rationalCurvePoint<nd, T>(u, deg_, knots_, cp_, w_, pt);
		}
		return pt;
	}

	void pointAndDerivatives(double u, std::vector<vecnt> &ptder) {
		if (!isRat_) {
			nurbs::curveDerivatives<nd, T>(u, deg_, knots_, cp_, 1, ptder);
		}
		else {
			//nurbs::curveDerivatives<nd, T>(u, deg_, knots_, cp_, w_, 1, ptder);
		}
	}

	vecnt tangent(double u) {
		std::vector<vecnt> tgt;
		if (!isRat_) {
			nurbs::curveDerivatives<nd, T>(u, deg_, knots_, cp_, 1, tgt);
		}
		else {
			// nurbs::rationalCurvePoint<nd, T>(u, deg_, knots_, cp_, w_, pt);
		}
		return tgt[1];
	}

	vecnt controlPoint(int i) const {
		return cp_[i];
	}

	void setControlPoint(int i, vecnt pt) {
		cp_[i] = pt;
	}

	std::vector<vecnt> controlPoints() const {
		return cp_;
	}

	int numControlPoints() const {
		return cp_.size();
	}

	void setRational(bool isRational) {
		isRat_ = isRational;
	}

	void setClampedAtStart(bool flag) {
		if (flag && !isClamp0_) {
			double startKnot = knots_[0];
			knots_.insert(knots_.begin(), deg_, startKnot);
		}
		if (!flag && isClamp0_) {
			knots_.erase(knots_.begin(), knots_.begin() + deg_);
		}
		isClamp0_ = flag;
	}

	void setClampedAtEnd(bool flag) {
		if (flag && !isClamp1_) {
			double endKnot = knots_[knots_.size() - 1];
			knots_.insert(knots_.end(), deg_, endKnot);
		}
		if (!flag && isClamp1_) {
			knots_.erase(knots_.end() - deg_, knots_.end());
		}
		isClamp1_ = flag;
	}

	void setClosed(bool flag) {
		if (flag && !isClosed_) {
			// wrap knots and duplicate last controlpoint
		}
		if (!flag && isClosed_) {
			// remove duplicate control point and unwrap knots
		}
	}

private:
	unsigned int deg_;
	std::vector<double> knots_;
	std::vector<vecnt> cp_;
	std::vector<T> w_;
	bool isRat_;
	bool isClamp0_, isClamp1_;
	bool isClosed_;
};

// Typedefs for ease of use
typedef NurbsCurve<2, float> NurbsCurve2f;
typedef NurbsCurve<2, double> NurbsCurve2d;
typedef NurbsCurve<3, float> NurbsCurve3f;
typedef NurbsCurve<3, double> NurbsCurve3d;

} // namespace ospl