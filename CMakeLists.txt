message(STATUS "=== PYNG CMakeLists LOADED: MARKER_ABC123 ===")
include_directories(SYSTEM /usr/include)
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

message(STATUS "NGSPICE include dirs: ${NGSPICE_INCLUDE_DIRS}")
message(STATUS "NGSPICE library dirs: ${NGSPICE_LIBRARY_DIRS}")
message(STATUS "NGSPICE libraries:    ${NGSPICE_LIBRARIES}")
message(STATUS "NGSPICE cflags:       ${NGSPICE_CFLAGS_OTHER}")

# ---- Core C++ library ------------------------------------------------------
add_library(pyng_core STATIC
        SpiceEngine.cpp
        montecarlo.cpp
)

# This is the key: explicitly add include dirs from pkg-config
target_include_directories(pyng_core
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        PRIVATE
        ${NGSPICE_INCLUDE_DIRS}
)

# Also explicitly add library search dirs + libs from pkg-config
if(NGSPICE_LIBRARY_DIRS)
    target_link_directories(pyng_core PRIVATE ${NGSPICE_LIBRARY_DIRS})
endif()

target_link_libraries(pyng_core PRIVATE ${NGSPICE_LIBRARIES})

# Some pkg-config packages require extra flags (rare here, but correct)
if(NGSPICE_CFLAGS_OTHER)
    target_compile_options(pyng_core PRIVATE ${NGSPICE_CFLAGS_OTHER})
endif()

# ---- Python extension module (import name: pyng) ---------------------------
pybind11_add_module(pyng MODULE
        python_external.cpp
)

target_link_libraries(pyng PRIVATE pyng_core)

# The Python module may also need ngspice on the link line depending on how
# your core library is built/linked on the platform
if(NGSPICE_LIBRARY_DIRS)
    target_link_directories(pyng PRIVATE ${NGSPICE_LIBRARY_DIRS})
endif()
target_link_libraries(pyng PRIVATE ${NGSPICE_LIBRARIES})

target_include_directories(pyng PRIVATE ${NGSPICE_INCLUDE_DIRS})
if(NGSPICE_CFLAGS_OTHER)
    target_compile_options(pyng PRIVATE ${NGSPICE_CFLAGS_OTHER})
endif()

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
