project(WaveOrthotope)
cmake_minimum_required(VERSION 3.19)

# Use C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS FALSE)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

add_executable(wavesolve_cmake src/main.cpp src/wave_orthotope.cpp )
