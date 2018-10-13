cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 14)

option(BUILD_TESTS "Build unit tests" ON)

project(tinynurbs LANGUAGES CXX)

find_package(glm 0.9.9 REQUIRED)

if (glm_FOUND)
    message("Found glm at ${glm_DIR}")
endif()

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
