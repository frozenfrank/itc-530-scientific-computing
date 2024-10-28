cmake_minimum_required(VERSION 3.20)
project(wavesolve CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS FALSE)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

add_executable(wavesolve_serial wavesolve_serial.cpp binary_io.hpp wave_orthotope.hpp)
