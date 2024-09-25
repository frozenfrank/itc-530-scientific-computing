cmake_minimum_required(VERSION 3.20)
project(wavesolve CXX)

# Use C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS FALSE)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Build wavesolve
add_executable(wavesolve wavesolve.cpp WaveOrthotope.hpp)