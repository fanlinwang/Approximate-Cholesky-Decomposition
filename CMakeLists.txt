project(Approximate_Cholesky_for_SDDM)

cmake_minimum_required(VERSION 2.8)

# C++11
set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -ffast-math")
find_package(Boost 1.45.0) 

# compiler options
add_definitions(-march=native)

add_subdirectory(src)
add_subdirectory(test)

include_directories( ${Boost_INCLUDE_DIRS} )