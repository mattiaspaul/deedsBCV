cmake_minimum_required(VERSION 2.8)
project(deedsBCV)

include_directories("./src/")

SET(CMAKE_CXX_FLAGS  "-std=c++11 -O3 -fopenmp -mavx2 -msse4.2")

add_executable(linear src/linearBCV.cpp)
target_link_libraries(linear -lz)

add_executable(deeds src/deedsBCV0.cpp)
target_link_libraries(deeds -lz)

add_executable(applyBCV src/applyBCV.cpp)
target_link_libraries(applyBCV -lz)
