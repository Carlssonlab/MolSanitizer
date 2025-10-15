cmake_minimum_required(VERSION 3.16)
project(msani_stereoisomers)

# Set C++ standard based on RDKit version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# For static linking and portability - MUST be set before any targets are created
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(BUILD_SHARED_LIBS OFF)

# Find Python and pybind11
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

# Find RDKit
set(RDKit_ROOT_DIR "$ENV{RDBASE}")
if(NOT RDKit_ROOT_DIR)
    set(RDKit_ROOT_DIR "$ENV{CONDA_PREFIX}")
endif()

message(STATUS "RDKit root directory: ${RDKit_ROOT_DIR}")

# Find RDKit include directory
find_path(RDKIT_INCLUDE_DIR GraphMol/ROMol.h
    PATHS ${RDKit_ROOT_DIR}/include/rdkit
          ${RDKit_ROOT_DIR}/include
          /usr/include/rdkit
          /usr/local/include/rdkit
    NO_DEFAULT_PATH
)

# Find RDKit libraries
set(RDKIT_LIB_NAMES
    RDKitGraphMol
    RDKitSmilesParse
    RDKitRDGeneral
    RDKitDistGeomHelpers
    RDKitForceFieldHelpers
    RDKitForceField
)

set(RDKIT_LIBRARIES "")
foreach(lib_name ${RDKIT_LIB_NAMES})
    find_library(${lib_name}_LIB ${lib_name}
        PATHS ${RDKit_ROOT_DIR}/lib
              /usr/lib
              /usr/local/lib
        NO_DEFAULT_PATH
    )
    if(${lib_name}_LIB)
        list(APPEND RDKIT_LIBRARIES ${${lib_name}_LIB})
        message(STATUS "Found ${lib_name}: ${${lib_name}_LIB}")
    else()
        message(WARNING "Could not find ${lib_name}")
    endif()
endforeach()

# Find Boost (required by RDKit)
find_package(Boost REQUIRED COMPONENTS system filesystem serialization)

# Create static library for stereoisomer enumeration
add_library(stereoisomer_lib STATIC
    EnumerateStereoisomers.cpp
    Flippers.cpp
)

target_include_directories(stereoisomer_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${RDKIT_INCLUDE_DIR}
)

target_link_libraries(stereoisomer_lib PUBLIC
    ${RDKIT_LIBRARIES}
    Boost::system
    Boost::filesystem
    Boost::serialization
)

# Create the Python module
pybind11_add_module(msani_stereoisomers
    pybind_interface.cpp
)

target_link_libraries(msani_stereoisomers PRIVATE
    stereoisomer_lib
)

target_compile_definitions(msani_stereoisomers PRIVATE VERSION_INFO="0.4.0")

# Install the Python module
install(TARGETS msani_stereoisomers
    LIBRARY DESTINATION .
    COMPONENT python
)
