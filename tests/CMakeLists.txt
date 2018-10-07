add_executable(unittests catch.hpp test_curve.cpp)
add_test(NAME unittests COMMAND unittests)
