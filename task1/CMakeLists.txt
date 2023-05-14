cmake_minimum_required(VERSION 3.25)
project(san9_parallel)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp")

add_executable(san9_parallel main.cpp)