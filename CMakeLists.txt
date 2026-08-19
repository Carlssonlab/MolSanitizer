cmake_minimum_required(VERSION 3.20)
project(MolSanitizer)

# Installed wheels must never inherit build-machine or Conda library paths.
# auditwheel will add relative wheel-local paths only when a distribution build
# intentionally vendors remaining shared dependencies.
set(CMAKE_SKIP_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

include(cmake/MSANINative.cmake)

# Find Python and pybind11 (shared across all modules)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Build conformer generation C++ extension
add_subdirectory(msani/conformers/cpp)

# Build stereoisomer enumeration C++ extension  
add_subdirectory(msani/moltransform/cpp)

# Build the vendored AMSOLcpp Python package.
add_subdirectory(msani/amsol)
