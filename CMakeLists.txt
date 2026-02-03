cmake_minimum_required(VERSION 3.22)
project(pyng LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(pyng
        SpiceEngine.cpp
        SpiceEngine.hpp
        main.cpp
        montecarlo.cpp
        montecarlo.hpp
)

find_path(NGSPICE_INCLUDE_DIR
        NAMES ngspice/sharedspice.h
)

find_library(NGSPICE_LIBRARY
        NAMES ngspice libngspice
)

if(NOT NGSPICE_INCLUDE_DIR OR NOT NGSPICE_LIBRARY)
    message(FATAL_ERROR "Could not find ngspice (sharedspice.h / libngspice). Set NGSPICE_INCLUDE_DIR and NGSPICE_LIBRARY.")
endif()

target_include_directories(pyng PRIVATE ${NGSPICE_INCLUDE_DIR})
target_link_libraries(pyng PRIVATE ${NGSPICE_LIBRARY})