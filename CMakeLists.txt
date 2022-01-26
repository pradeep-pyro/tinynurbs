cmake_minimum_required(VERSION 3.5)
project(tinynurbs LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

set(BUILD_TESTS ON CACHE BOOL "Build unit tests")

find_package(GLM 0.9.9 REQUIRED MODULE)


set(HEADER_FILES 
    include/tinynurbs/tinynurbs.h
    include/tinynurbs/core/basis.h
    include/tinynurbs/core/check.h
    include/tinynurbs/core/curve.h
    include/tinynurbs/core/evaluate.h
    include/tinynurbs/core/modify.h
    include/tinynurbs/core/surface.h
    include/tinynurbs/io/obj.h
    include/tinynurbs/util/util.h
    include/tinynurbs/util/array2.h
)
source_group("Header Files" FILES ${HEADER_FILES})
source_group("CMake Files" FILES CMakeLists.txt)

add_library(tinynurbs INTERFACE)
add_library(tinynurbs::tinynurbs ALIAS tinynurbs )
target_link_libraries(tinynurbs INTERFACE glm)

target_include_directories(tinynurbs INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

add_custom_target(tinynurbs_dummy SOURCES ${HEADER_FILES} CMakeLists.txt)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
