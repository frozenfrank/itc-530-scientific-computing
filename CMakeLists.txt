project(WaveOrthotope)
cmake_minimum_required(VERSION 3.19)

# Use C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS FALSE)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

set(COMMON_INCLUDES src/wave_orthotope.hpp src/wave_orthotope.cpp src/binary_io.hpp)

# Executables
add_executable(wavesolve src/main.cpp src/wave_orthotope.hpp src/wave_orthotope.cpp)
add_executable(wavesolve_serial src/main_serial.cpp ${COMMON_INCLUDES})
