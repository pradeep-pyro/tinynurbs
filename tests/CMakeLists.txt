add_executable(unittests catch.hpp test_main.cpp
                         test_curve.cpp test_rational_curve.cpp
                         test_surface.cpp test_rational_surface.cpp)
add_test(NAME unittests COMMAND unittests)
