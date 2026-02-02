cmake_minimum_required(VERSION 4.1)
project(pyng)

set(CMAKE_CXX_STANDARD 23)

include_directories(.)

add_executable(pyng
        engine.cpp
        engine.hpp
        main.cpp
        montecarlo.cpp
        montecarlo.hpp)

target_link_libraries(pyng ngspice)