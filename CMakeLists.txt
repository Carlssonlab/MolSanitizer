cmake_minimum_required(VERSION 3.20)
project(MolSanitizer)

# Find Python and pybind11 (shared across all modules)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Build conformer generation C++ extension
add_subdirectory(msani/conformers/cpp)

# Build stereoisomer enumeration C++ extension  
add_subdirectory(msani/moltransform/cpp)

# Build the vendored AMSOLcpp Python package.
add_subdirectory(msani/amsol)
