cmake_minimum_required(VERSION 3.8)
project(math_lib VERSION 1.0.0 LANGUAGES CXX)

set(LIBRARY_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)

add_subdirectory(lib)

enable_testing()
add_subdirectory(tests)