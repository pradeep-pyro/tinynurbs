cmake_minimum_required(VERSION 3.1.0)
set(CMAKE_CXX_STANDARD 11)

project(nurbstk)

add_definitions(-DGLM_DEFAULT_CTOR)

message("${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
find_package(glm REQUIRED)
message("Using glm from: ${glm_INCLUDE_DIR}")

file(GLOB INCLUDES "nurbstk/*.h")
file(GLOB SOURCES "nurbstk/*.cpp")

# add_library(nurbstk STATIC ${SOURCES})
add_executable(nurbstk main.cpp ${INCLUDES} ${SOURCES})
target_include_directories(nurbstk PUBLIC ${glm_INCLUDE_DIR})

