cmake_minimum_required(VERSION 3.20)
# Match official CPython's DLL runtime. The Windows RDKit/Boost SDK must also
# be built with /MD (static libraries do not imply a static C runtime).
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
project(MolSanitizer)

# Installed wheels must never inherit build-machine or Conda library paths.
# auditwheel will add relative wheel-local paths only when a distribution build
# intentionally vendors remaining shared dependencies.
set(CMAKE_SKIP_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

include(cmake/MSANINative.cmake)

# Find Python and pybind11 (shared across all modules)
# Extensions need module headers/targets, not a Python embedding library.
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)

# pip installs pybind11's CMake config under site-packages, which is not a
# standard CMake prefix. Query the selected build interpreter, not the SDK.
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -m pybind11 --cmakedir
    RESULT_VARIABLE _msani_pybind11_status
    OUTPUT_VARIABLE _msani_pybind11_cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(_msani_pybind11_status EQUAL 0 AND
   IS_DIRECTORY "${_msani_pybind11_cmake_dir}")
    list(APPEND CMAKE_PREFIX_PATH "${_msani_pybind11_cmake_dir}")
endif()
find_package(pybind11 CONFIG REQUIRED)

# Native implementation and Python bindings share one source entry point.
add_subdirectory(msani/cpp)
