# tinynurbs

This is a lightweight C++ library for Non-Uniform Rational B-Spline curves and surfaces. The API is simple to use and the code is quite readable while being efficient.

Current features mainly focus on point & derivative evaluations, and Wavefront OBJ file reading/writing (for single NURBS surface). The library is under active development.

## Dependencies

* [glm] (version 0.9.9 where PR #584 is merged is required since tinynurbs uses the `glm::vec<dim, T>` type)

[glm]: https://github.com/g-truc/glm