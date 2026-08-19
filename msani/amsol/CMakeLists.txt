cmake_minimum_required(VERSION 3.20)

project(AMSOLcpp VERSION 0.1.0 LANGUAGES CXX)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

add_library(amsolcpp_core STATIC
    src/api.cpp
    src/coordinates.cpp
    src/elements.cpp
    src/integrals.cpp
    src/matrix.cpp
    src/parser.cpp
    src/scf.cpp
    src/solvation.cpp
    src/writer.cpp
)
msani_configure_native_static_target(amsolcpp_core)
target_compile_features(amsolcpp_core PUBLIC cxx_std_20)
target_include_directories(amsolcpp_core
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
)
if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    # Preserve AMSOLcpp's audited floating-point behavior across builds.
    target_compile_options(amsolcpp_core PRIVATE
        -fno-fast-math
        -ffp-contract=off
    )
endif()

pybind11_add_module(_amsolcpp MODULE python/bindings.cpp)
msani_configure_native_extension(_amsolcpp)
target_compile_features(_amsolcpp PRIVATE cxx_std_20)
target_link_libraries(_amsolcpp PRIVATE amsolcpp_core)

# Enable the direct RDKit Mol API when RDKit's development files are present.
# Conda normally provides a CMake package; the manual lookup also supports
# environments that ship only the headers and shared libraries.
find_package(RDKit CONFIG QUIET)
if(NOT RDKit_FOUND)
    find_path(RDKIT_INCLUDE_DIR GraphMol/ROMol.h
        PATHS
            $ENV{RDBASE}/include/rdkit
            $ENV{CONDA_PREFIX}/include/rdkit
            /usr/include/rdkit
            /usr/local/include/rdkit
    )
endif()

if(RDKit_FOUND)
    msani_validate_rdkit_linkage(
        RDKit::GraphMol RDKit::RDGeometryLib RDKit::RDGeneral
    )
    target_compile_definitions(_amsolcpp PRIVATE AMSOLCPP_WITH_RDKIT)
    target_link_libraries(_amsolcpp PRIVATE
        RDKit::GraphMol
        RDKit::RDGeometryLib
        RDKit::RDGeneral
    )
elseif(RDKIT_INCLUDE_DIR)
    target_compile_definitions(_amsolcpp PRIVATE AMSOLCPP_WITH_RDKIT)
    target_include_directories(_amsolcpp PRIVATE ${RDKIT_INCLUDE_DIR})
    foreach(_rdkit_library RDKitGraphMol RDKitRDGeneral RDKitRDGeometryLib)
        find_library(${_rdkit_library}_LIB NAMES ${_rdkit_library}
            PATHS
                $ENV{RDBASE}/lib
                $ENV{CONDA_PREFIX}/lib
                /usr/lib
                /usr/local/lib
        )
        if(${_rdkit_library}_LIB)
            target_link_libraries(_amsolcpp PRIVATE ${${_rdkit_library}_LIB})
        endif()
    endforeach()
else()
    message(STATUS
        "RDKit development files not found; amsolcpp direct RDKit functions "
        "will report that RDKit support is unavailable")
endif()

install(TARGETS _amsolcpp
    LIBRARY DESTINATION amsolcpp COMPONENT python
    RUNTIME DESTINATION amsolcpp COMPONENT python
)
install(FILES python/amsolcpp/__init__.py
    DESTINATION amsolcpp
    COMPONENT python
)
