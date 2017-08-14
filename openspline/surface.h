/*
@file openspline/surface.h
@author Pradeep Kumar Jayaraman <pradeep.pyro@gmail.com>

The Surface class represents an abstract 3D paramteric surface.

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE.txt file.
*/

#pragma once
#include "glm/glm.hpp"

namespace ospl {
template <int nd, typename T>
class Surface {
	typedef glm::vec<nd, T> vecnt;
public:
	virtual unsigned int degreeU() const = 0;
	virtual unsigned int degreeV() const = 0;
	virtual vecnt point(double u, double v) const = 0;
	virtual void derivatives(double u, double v, int nDers,
		std::vector<std::vector<vecnt>> &ders) const = 0;
	virtual void tangent(double u, double v,
		vecnt &du, vecnt &dv, bool normalize=true) const = 0;
	virtual vecnt normal(double u, double v, bool normalize=true) const = 0;
};

typedef Surface<3, double> Surface3d;
typedef Surface<3, float> Surface3f;

} // namespace ospl