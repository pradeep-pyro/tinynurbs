# tinynurbs

This is a lightweight header-only C++ library for Non-Uniform Rational B-Spline curves and surfaces. The API is simple to use and the code is readable while being efficient.

Some of the main features include:

- Point & derivative evaluations upto any order
- Support for clamped and unclamped curves/surfaces
- Knot insertion
- Wavefront OBJ file I/O

The library is under active development.

## Usage

The entire API consists of free functions named as `curve*` and `surface*` for curve- and surface-based  functionalities.

`Examples coming soon...`

## Dependencies

* [glm] (version 0.9.9 where PR #584 is merged is required since tinynurbs uses the `glm::vec<dim, T>` type)

[glm]: https://github.com/g-truc/glm