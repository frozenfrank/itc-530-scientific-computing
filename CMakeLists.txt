cmake_minimum_required(VERSION 3.20)
project(wavesolve CXX)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS FALSE)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG" CACHE STRING "Flags used by the CXX compiler during RELEASE builds")

add_executable(wavesolve_serial wavesolve_serial.cpp binary_io.hpp)
add_executable(wavesolve_openmp wavesolve_openmp.cpp binary_io.hpp)
target_link_libraries(wavesolve_openmp PRIVATE OpenMP::OpenMP_CXX)
add_executable(optimize optimize.cpp)
target_link_libraries(optimize PRIVATE OpenMP::OpenMP_CXX)
