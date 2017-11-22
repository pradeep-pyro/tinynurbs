#pragma once

#include <vector>
#include <array>

namespace nurbstk {

/**
// A simple class for representing 2D runtime arrays.
*/
template <typename T>
class array2 {
public:
	array2(size_t nRows, size_t nCols, T fillValue = 0.0)
		: rows(nRows), cols(nCols) {
		data.resize(rows * cols, fillValue);
	}
	T operator()(size_t row, size_t col) const {
		return data[row*cols + col];
	}
	T& operator()(size_t row, size_t col) {
		return data[row*cols + col];
	}
private:
	size_t rows, cols;
	std::vector<T> data;
};

} // namespace nurbstk