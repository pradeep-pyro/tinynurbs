#pragma once

#include <glm/glm.hpp>

#include <vector>

namespace nurbstk {
namespace util {

/**
// A simple class for representing 2D runtime arrays.
*/
template <typename T>
class array2 {
public:
    array2() {
    }
    array2(size_t rows, size_t cols, T default_value = T()) {
        resize(rows, cols, default_value);
    }
    void resize(size_t rows, size_t cols, T val=T()) {
        data_.resize(rows * cols, val);
        rows_ = rows;
        cols_ = cols;
    }
    T operator()(size_t row, size_t col) const {
        return data_[row*cols_ + col];
    }
    T& operator()(size_t row, size_t col) {
        return data_[row*cols_ + col];
    }
    T operator[](size_t idx) const {
        return data_[idx];
    }
    T& operator[](size_t idx) {
        return data_[idx];
    }
    size_t rows() const {
        return rows_;
    }
    size_t cols() const {
        return cols_;
    }
    size_t size() const {
        return data_.size();
    }
    T * data() {
        return &data_[0];
    }
private:
    size_t rows_, cols_;
    std::vector<T> data_;
};

/**
Convert an nd point in homogenous coordinates to an (n-1)d point in cartesian
coordinates by perspective division
@param[in] pt Point in homogenous coordinates
@return Input point in cartesian coordinates
*/
template<int nd, typename T>
inline glm::vec<nd - 1, T> homogenousToCartesian(const glm::vec<nd, T> &pt) {
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
inline glm::vec<nd + 1, T> cartesianToHomogenous(const glm::vec<nd, T> &pt, T w) {
    return glm::vec<nd + 1, T>(pt * w, w);
}

/**
Convert an (n+1)d point to an nd point without perspective division
by truncating the last dimension
@param[in] pt Point in homogenous coordinates
@return Input point in cartesian coordinates
*/
template<int nd, typename T>
inline glm::vec<nd - 1, T> truncateHomogenous(const glm::vec<nd, T> &pt) {
    return glm::vec<nd - 1, T>(pt);
}

/**
Compute the binomial coefficient (nCk) using the formula
\product_{i=0}^k (n + 1 - i) / i
*/
unsigned int binomial(unsigned int n, unsigned int k);

template <typename T>
inline bool close(T a, T b, double eps = std::numeric_limits<T>::epsilon()) {
    return (std::abs(a - b) < eps) ? true : false;
}


} // namespace util
} // namespace nurbstk
