cmake_minimum_required(VERSION 3.16)
project(fem_advection)

set(CMAKE_CXX_STANDARD 17)

add_executable(fem_advection main.cpp solve.cpp solve.h)