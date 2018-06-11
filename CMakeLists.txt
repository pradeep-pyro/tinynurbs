cmake_minimum_required(VERSION 3.6)
set(CMAKE_CXX_STANDARD 11)

project(nurbstk)

set(glm_dir "/Users/pradeep/libs/glm" CACHE STRING "Path to the glm library" FORCE)
message("Using glm from: ${glm_dir}")
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
target_include_directories(nurbstk PUBLIC ${glm_dir})
