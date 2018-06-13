cmake_minimum_required(VERSION 3.1.0)
set(CMAKE_CXX_STANDARD 11)

project(tinynurbs)

message("${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
find_package(glm REQUIRED)

file(GLOB_RECURSE INCLUDES "tinynurbs/*.h")
file(GLOB_RECURSE SOURCES "tinynurbs/*.cpp")

add_executable(${PROJECT_NAME} main.cpp ${INCLUDES} ${SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE ${glm_INCLUDE_DIR})
