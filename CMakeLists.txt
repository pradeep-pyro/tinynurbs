cmake_minimum_required(VERSION 3.1.0)
set(CMAKE_CXX_STANDARD 11)

project(nurbstk)

message("${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
find_package(glm REQUIRED)
message("Using glm from: ${glm_INCLUDE_DIR}")

set(SOURCES nurbstk/basis.h
            nurbstk/evaluate.h
            nurbstk/evaluate.cpp
            nurbstk/knots.h
            nurbstk/curve.h nurbstk/curve.cpp
            nurbstk/surface.h nurbstk/surface.cpp
            nurbstk/util.h nurbstk/util.cpp
            nurbstk/refine.h
            nurbstk/io.h)

# add_library(nurbstk STATIC ${SOURCES})
add_executable(nurbstk main.cpp ${SOURCES})
target_include_directories(nurbstk PUBLIC ${glm_INCLUDE_DIR})
