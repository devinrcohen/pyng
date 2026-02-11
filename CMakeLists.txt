cmake_minimum_required(VERSION 3.22)
project(pyng LANGUAGES C CXX)

option(PYNG_BUILD_CLI "Build CLI executable (main.cpp)" ON)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ---- Python + pybind11 -----------------------------------------------------
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# ---- ngspice via pkg-config ------------------------------------------------
find_package(PkgConfig REQUIRED)
pkg_check_modules(NGSPICE REQUIRED ngspice)

message(STATUS "ngspice via pkg-config: version ${NGSPICE_VERSION}")
message(STATUS "NGSPICE include dirs: ${NGSPICE_INCLUDE_DIRS}")
message(STATUS "NGSPICE library dirs: ${NGSPICE_LIBRARY_DIRS}")
message(STATUS "NGSPICE libs: ${NGSPICE_LIBRARIES}")

# Resolve the actual library file (preferred) so we don't rely on link search paths.
find_library(NGSPICE_LIBRARY_FILE
        NAMES ngspice libngspice
        HINTS ${NGSPICE_LIBRARY_DIRS}
        PATHS /usr/lib64 /usr/lib /lib64 /lib /usr/local/lib
)
if(NOT NGSPICE_LIBRARY_FILE)
    message(FATAL_ERROR "Could not resolve libngspice via find_library(). "
            "pkg-config said libs='${NGSPICE_LIBRARIES}' dirs='${NGSPICE_LIBRARY_DIRS}'")
endif()
message(STATUS "NGSPICE library file: ${NGSPICE_LIBRARY_FILE}")

# ---- Core C++ library ------------------------------------------------------
add_library(pyng_core STATIC
        SpiceEngine.cpp
        montecarlo.cpp
)

target_include_directories(pyng_core
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        # Still add include dirs normally (useful for IDEs)
        ${NGSPICE_INCLUDE_DIRS}
)

# Force include dirs even if CMake thinks they're "implicit".
# This is the key fix for conda toolchain compilers that don't search /usr/include.
foreach(_inc IN LISTS NGSPICE_INCLUDE_DIRS)
    target_compile_options(pyng_core PRIVATE "-I${_inc}")
endforeach()

target_link_libraries(pyng_core
        PUBLIC
        ${NGSPICE_LIBRARY_FILE}
)

# ---- Python extension module (pyng) ---------------------------------------
pybind11_add_module(pyng MODULE
        python_external.cpp
)

target_link_libraries(pyng PRIVATE pyng_core)

# ---- Optional CLI ----------------------------------------------------------
if(PYNG_BUILD_CLI)
    add_executable(pyng_cli main.cpp)
    target_link_libraries(pyng_cli PRIVATE pyng_core)
endif()

# ---- Warnings --------------------------------------------------------------
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    target_compile_options(pyng_core PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(pyng      PRIVATE -Wall -Wextra -Wpedantic)
endif()
