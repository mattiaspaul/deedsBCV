cmake_minimum_required(VERSION 3.19)
project(deedsBCV)

find_package ( OpenMP REQUIRED )

include_directories("./src/")

add_executable(linear src/linearBCV.cpp)
target_link_libraries(linear PRIVATE OpenMP::OpenMP_C -lz)

add_executable(deeds src/deedsBCV0.cpp)
target_link_libraries(deeds PRIVATE OpenMP::OpenMP_C -lz)

add_executable(applyBCV src/applyBCV.cpp)
target_link_libraries(applyBCV PRIVATE OpenMP::OpenMP_C -lz)
