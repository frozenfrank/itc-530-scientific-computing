cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -fopenmp" CACHE STRING "Flags used by the CXX compiler during RELEASE builds")
project(wavesolve CXX)
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)

# Use C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS FALSE)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Build wavesolve
add_executable(wavesolve_openmp wavesolve_openmp.cpp  WaveOrthotope.hpp )
target_link_libraries(wavesolve_openmp PRIVATE OpenMP::OpenMP_CXX)
add_executable(optimize optimize.cpp )
target_link_libraries(optimize PRIVATE OpenMP::OpenMP_CXX)
add_executable(wavesolve_thread wavesolve_thread.cpp WaveOrthotope.hpp)
target_link_libraries(wavesolve_thread PRIVATE Threads::Threads)