cmake_minimum_required(VERSION 3.24)
project(task2)

set(CMAKE_CXX_STANDARD 17)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_executable(task2 main.cpp)

target_link_libraries(task2 ${MPI_CXX_LIBRARIES})