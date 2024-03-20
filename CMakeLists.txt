cmake_minimum_required(VERSION 3.4)
project(polynomial)
add_subdirectory(pybind11)
pybind11_add_module(polynomial main.cpp)