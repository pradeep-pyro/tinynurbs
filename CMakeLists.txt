cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 14)

project(tinynurbs LANGUAGES CXX)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
