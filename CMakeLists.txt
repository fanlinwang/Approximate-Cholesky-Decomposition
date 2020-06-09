project(Approximate_Cholesky_for_SDDM)

cmake_minimum_required(VERSION 2.8)

# C++11
set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -ffast-math -fno-tree-vectorize -pg -g")
find_package(Boost 1.45.0) 
set(CMAKE_BUILD_TYPE Debug)

# compiler options
add_definitions(-march=native)

add_subdirectory(src)
add_subdirectory(test)

include_directories( ${Boost_INCLUDE_DIRS} )