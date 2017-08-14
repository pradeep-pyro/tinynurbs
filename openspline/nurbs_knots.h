/*
@file openspline/nurbs_knots.h
@author Pradeep Kumar Jayaraman <pradeep.pyro@gmail.com>

Helper functions for creating and modifying knot vectors.

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE.txt file.
*/

#pragma once

#include <vector>

bool makeUniformKnotVector(unsigned int degree, size_t nCtrlPts,
	std::vector<double> &knots);

bool makeClampedUniformKnotVector(unsigned int degree, size_t nCtrlPts,
	std::vector<double> &knots);

void clampKnotVector(unsigned int degree, std::vector<double> &knots);

void clampKnotVectorLeft(unsigned int degree, std::vector<double> &knots);

void clampKnotVectorRight(unsigned int degree, std::vector<double> &knots);

bool isKnotVectorMonotonic(const std::vector<double> &knots);

bool isKnotVectorClosed(unsigned int degree, const std::vector<double> &knots);
