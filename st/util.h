#pragma once

#include <glm/glm.hpp>

#include <vector>

namespace nurbstk {
namespace util {

/**
Convert an nd point in homogenous coordinates to an (n-1)d point in cartesian
coordinates by perspective division
@param[in] pt Point in homogenous coordinates
@return Input point in cartesian coordinates
*/
template<int nd, typename T>
inline glm::vec<nd - 1, T> homogenousToCartesian(glm::vec<nd, T> pt) {
	return glm::vec<nd - 1, T>(pt / pt[pt.length() - 1]);
}

/**
Convert an nd point in cartesian coordinates to an (n+1)d point in homogenous
coordinates
@param[in] pt Point in cartesian coordinates
@param[in] w Weight
@return Input point in homogenous coordinates
*/
template<int nd, typename T>
inline glm::vec<nd + 1, T> cartesianToHomogenous(glm::vec<nd, T> pt, T w) {
	return glm::vec<nd + 1, T>(pt * w, w);
}

/**
Convert an (n+1)d point to an nd point without perspective division
by truncating the last dimension
@param[in] pt Point in homogenous coordinates
@return Input point in cartesian coordinates
*/
template<int nd, typename T>
inline glm::vec<nd - 1, T> truncateHomogenous(glm::vec<nd, T> pt) {
	return glm::vec<nd - 1, T>(pt);
}

/**
Compute the binomial coefficient (nCk) using the formula
\product_{i=0}^k (n + 1 - i) / i
*/
unsigned int binomial(unsigned int n, unsigned int k);

} // namespace util
} // namespace nurbstk