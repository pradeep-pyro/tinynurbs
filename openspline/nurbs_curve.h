/*
@file openspline/nurbs_curve.h
@author Pradeep Kumar Jayaraman <pradeep.pyro@gmail.com>

The NurbsCurve class represents an n-dimensional non-uniform rational curve.

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
		const std::vector<vecnt> &controlPoints) throw (std::logic_error) {
		if (knots.size() != controlPoints.size() + degree + 1) {
			throw std::logic_error("Mismatch between degree, knots and control \
		points");
			return;
		}
		deg_ = degree;
		knots_ = knots;
		cp_ = controlPoints;
		isRational_ = false;
	}

	NurbsCurve(unsigned int degree, const std::vector<double> &knots, 
		const std::vector<vecnt> &controlPoints, const std::vector<T> &weights)
		throw (std::logic_error) {
		if (knots.size() != controlPoints.size() + degree + 1) {
			throw std::logic_error("Mismatch between degree, knots and control \
				points");
			return;
		}
		deg_ = degree;
		knots_ = knots;
		cp_ = controlPoints;
		w_ = weights;
		isRational_ = true;
	}
	
	vecnt pointAt(T u) {
		vecnt pt;
		if (!isRational_) {
			nurbs::curvePoint(u, deg_, knots_, cp_, pt);
		}
		else {
			nurbs::rationalCurvePoint(u, deg_, knots_, cp_, w_, pt);
		}
		return pt;
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
		isRational_ = isRational;
	}

	void setClampedAtStart(bool flag) {
		if (flag && !isClamped0_) {
			double startKnot = knots_[0];
			knots_.insert(knots_.begin(), deg_, startKnot);
		}
		if (!flag && isClamped0_) {
			knots_.erase(knots.begin(), knots.begin() + deg_);
		}
		isClamped0_ = flag;
	}

	void setClampedAtEnd(bool flag) {
		if (flag && !isClamped1_) {
			double endKnot = knots_[knots_.size() - 1];
			knots_.insert(knots_.end(), deg_, endKnot);
		}
		if (!flag && isClamped1_) {
			knots_.erase(knots.end() - deg_, knots.end());
		}
		isClamped1_ = flag;
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
	bool isRational_;
	bool isClamped0_, isClamped1_;
	bool isClosed_;
};


} // namespace ospl