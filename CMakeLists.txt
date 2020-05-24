project(Approximate_Cholesky_for_SDDM)

cmake_minimum_required(VERSION 2.8)

# C++11
set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -ffast-math")

# compiler options
add_definitions(-march=native)

add_subdirectory(src)
add_subdirectory(test)

set(Boost_USE_STATIC_LIBS        ON)  # only find static libs
set(Boost_USE_DEBUG_LIBS         OFF) # ignore debug libs and
set(Boost_USE_RELEASE_LIBS       ON)  # only find release libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

find_package(Boost)
IF (Boost_FOUND)
    message(STATUS "${Boost_INCLUDE_DIRS}")
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

ENDIF()