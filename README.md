# tinynurbs

This is a lightweight header-only C++14 library for Non-Uniform Rational B-Spline curves and surfaces. The API is simple to use and the code is readable while being efficient.

Some of the main features include:

- Supports non-rational and rational curves and surfaces of any order
- Supports planar and space curves
- Point & derivative evaluations upto any order
- Knot insertion, splitting without affecting the shape
- Wavefront OBJ format I/O

The library is under active development.

## Dependencies

- [glm] (version 0.9.9 where PR #584 is merged is required since tinynurbs uses the `glm::vec<dim, T>` type)
- C++14 compliant compiler

## Usage

The entire API consists of free functions named as `curve*` and `surface*` for curve- and surface-based  functionalities.
Examples and documentation can be found here:


## References

Some useful reference material and NURBS implementations:

-
-

[glm]: https://github.com/g-truc/glm
