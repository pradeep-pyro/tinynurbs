add_executable(unittests catch.hpp test_main.cpp
                         test_curve.cpp test_rational_curve.cpp
                         test_surface.cpp test_rational_surface.cpp)
target_include_directories(unittests PRIVATE ${GLM_INCLUDE_DIRS})
target_link_libraries(unittests PUBLIC tinynurbs::tinynurbs)
add_test(NAME unittests COMMAND unittests)
