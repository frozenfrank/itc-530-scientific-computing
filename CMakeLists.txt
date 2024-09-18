project(WaveOrthotope)
cmake_minimum_required(VERSION 3.19)

set(CMAKE_CXX_FLAGS "-std=c++20")
add_executable(wavesolve_cmake main.cpp wave_orthotope.cpp )
