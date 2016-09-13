cmake_minimum_required(VERSION 3.3)
project(BiQuad)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp BiQuad.cpp BiQuad.h)
add_executable(BiQuad ${SOURCE_FILES})