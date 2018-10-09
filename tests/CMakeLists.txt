add_executable(unittests catch.hpp test_curve.cpp test_surface.cpp)
add_test(NAME unittests COMMAND unittests)
