cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

set(BUILD_TESTS true CACHE BOOL "Build unit tests")
set(GLM_ROOT_DIR "" CACHE STRING "Root directory of GLM (>=0.9.9)")
message("Variable from cache: ${GLM_ROOT_DIR}")

project(tinynurbs LANGUAGES CXX)

find_package(GLM 0.9.9 REQUIRED MODULE)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
