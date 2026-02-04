cmake_minimum_required(VERSION 3.22)
project(pyng LANGUAGES C CXX)

option(PYNG_BUILD_CLI "Build the standalone CLI executable (main.cpp)" ON)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ---- Python + pybind11 (from conda or system Python) -----------------------
# Recommended when using conda:
#   -DPython3_EXECUTABLE=$CONDA_PREFIX/bin/python
#   -DPython3_ROOT_DIR=$CONDA_PREFIX
#   -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# ---- ngspice (sharedspice) ------------------------------------------------
set(_ngspice_hints)
if(DEFINED ENV{CONDA_PREFIX})
    list(APPEND _ngspice_hints
            "$ENV{CONDA_PREFIX}"
            "$ENV{CONDA_PREFIX}/include"
            "$ENV{CONDA_PREFIX}/lib"
    )
endif()
list(APPEND _ngspice_hints
        "/usr/local"
        "/usr"
        "/opt/homebrew"
        "/opt/local"
)

find_path(NGSPICE_INCLUDE_DIR
        NAMES ngspice/sharedspice.h
        HINTS ${_ngspice_hints}
        PATH_SUFFIXES include
)

find_library(NGSPICE_LIBRARY
        NAMES ngspice libngspice
        HINTS ${_ngspice_hints}
        PATH_SUFFIXES lib lib64
)

if(NOT NGSPICE_INCLUDE_DIR OR NOT NGSPICE_LIBRARY)
    message(FATAL_ERROR
            "Could not find ngspice (ngspice/sharedspice.h and libngspice).\n"
            "Set -DNGSPICE_INCLUDE_DIR=/path/to/include and -DNGSPICE_LIBRARY=/path/to/libngspice.\n"
            "Detected: NGSPICE_INCLUDE_DIR='${NGSPICE_INCLUDE_DIR}', NGSPICE_LIBRARY='${NGSPICE_LIBRARY}'"
    )
endif()

# ---- Core library (no Python) ---------------------------------------------
add_library(pyng_core STATIC
        SpiceEngine.cpp
        montecarlo.cpp
)

target_include_directories(pyng_core
        PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
        PRIVATE
        "${NGSPICE_INCLUDE_DIR}"
)

target_link_libraries(pyng_core
        PRIVATE
        "${NGSPICE_LIBRARY}"
)

# ---- Python extension module (name: pyng) ---------------------------------
# Produces pyng.<abi>.so so Python can `import pyng`
pybind11_add_module(pyng MODULE
        python_external.cpp
)

if(PYNG_BUILD_CLI)
    add_executable(pyng_cli main.cpp)
    target_link_libraries(pyng_cli PRIVATE pyng_core)
endif()

target_link_libraries(pyng
        PRIVATE
        pyng_core
        "${NGSPICE_LIBRARY}"
)

target_include_directories(pyng
        PRIVATE
        "${NGSPICE_INCLUDE_DIR}"
)

# Warnings (optional but helpful)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    target_compile_options(pyng_core PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(pyng      PRIVATE -Wall -Wextra -Wpedantic)
endif()

# ---- Optional CLI ----------------------------------------------------------
#if(PYNG_BUILD_CLI)
#    add_executable(pyng_cli main.cpp)
#    target_link_libraries(pyng_cli PRIVATE pyng_core)
#endif()
