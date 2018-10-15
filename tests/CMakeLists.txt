message("Variable from cache: ${GLM_ROOT_DIR}")
add_executable(unittests catch.hpp test_main.cpp
                         test_curve.cpp test_rational_curve.cpp
                         test_surface.cpp test_rational_surface.cpp)
target_include_directories(unittests PRIVATE ${GLM_INCLUDE_DIRS})
add_test(NAME unittests COMMAND unittests)
