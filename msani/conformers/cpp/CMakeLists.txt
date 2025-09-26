cmake_minimum_required(VERSION 3.16)
project(stochastic_sampling_combined)

# Find RDKit first to determine version for C++ standard
set(RDKit_ROOT_DIR "$ENV{RDBASE}")
if(NOT RDKit_ROOT_DIR)
    set(RDKit_ROOT_DIR "$ENV{CONDA_PREFIX}")
endif()

# Try to find RDKit version to determine C++ standard
find_path(RDKIT_INCLUDE_DIR GraphMol/ROMol.h
    PATHS ${RDKit_ROOT_DIR}/include/rdkit
          ${RDKit_ROOT_DIR}/include
          /usr/include/rdkit
          /usr/local/include/rdkit
    NO_DEFAULT_PATH)

if(RDKIT_INCLUDE_DIR)
    # Try to detect RDKit version to set appropriate C++ standard
    find_file(RDKIT_VERSION_HEADER
        NAMES Version.h version.h
        PATHS ${RDKIT_INCLUDE_DIR}/GraphMol
              ${RDKIT_INCLUDE_DIR}/RDKit
              ${RDKIT_INCLUDE_DIR}
        NO_DEFAULT_PATH)
    
    if(RDKIT_VERSION_HEADER)
        # Read the version header to determine RDKit version
        file(READ ${RDKIT_VERSION_HEADER} RDKIT_VERSION_CONTENT LIMIT 2048)
        
        # Try to extract version information
        if(RDKIT_VERSION_CONTENT MATCHES "RDKIT_VERSION[ \t]+\"([0-9]+)\\.([0-9]+)\\.([0-9]+)\"")
            set(RDKIT_VERSION_MAJOR ${CMAKE_MATCH_1})
            set(RDKIT_VERSION_MINOR ${CMAKE_MATCH_2})
            set(RDKIT_VERSION_PATCH ${CMAKE_MATCH_3})
            set(RDKIT_VERSION "${RDKIT_VERSION_MAJOR}.${RDKIT_VERSION_MINOR}.${RDKIT_VERSION_PATCH}")
            message(STATUS "Detected RDKit version: ${RDKIT_VERSION}")
            
            # Set C++ standard based on RDKit version
            # RDKit <= 2025.3.1: C++17, > 2025.3.1: C++20
            if(RDKIT_VERSION VERSION_LESS_EQUAL "2025.3.1")
                set(CMAKE_CXX_STANDARD 17)
                message(STATUS "Using C++17 for RDKit version ${RDKIT_VERSION}")
            else()
                set(CMAKE_CXX_STANDARD 20)
                message(STATUS "Using C++20 for RDKit version ${RDKIT_VERSION}")
            endif()
        elseif(RDKIT_VERSION_CONTENT MATCHES "RDKit_VERSION_MAJOR[ \t]+([0-9]+)")
            # Alternative version detection
            string(REGEX MATCH "RDKit_VERSION_MAJOR[ \t]+([0-9]+)" _ ${RDKIT_VERSION_CONTENT})
            set(RDKIT_VERSION_MAJOR ${CMAKE_MATCH_1})
            string(REGEX MATCH "RDKit_VERSION_MINOR[ \t]+([0-9]+)" _ ${RDKIT_VERSION_CONTENT})
            set(RDKIT_VERSION_MINOR ${CMAKE_MATCH_1})
            string(REGEX MATCH "RDKit_VERSION_PATCH[ \t]+([0-9]+)" _ ${RDKIT_VERSION_CONTENT})
            set(RDKIT_VERSION_PATCH ${CMAKE_MATCH_1})
            
            if(RDKIT_VERSION_MAJOR AND RDKIT_VERSION_MINOR AND RDKIT_VERSION_PATCH)
                set(RDKIT_VERSION "${RDKIT_VERSION_MAJOR}.${RDKIT_VERSION_MINOR}.${RDKIT_VERSION_PATCH}")
                message(STATUS "Detected RDKit version: ${RDKIT_VERSION}")
                
                if(RDKIT_VERSION VERSION_LESS_EQUAL "2025.3.1")
                    set(CMAKE_CXX_STANDARD 17)
                    message(STATUS "Using C++17 for RDKit version ${RDKIT_VERSION}")
                else()
                    set(CMAKE_CXX_STANDARD 20)
                    message(STATUS "Using C++20 for RDKit version ${RDKIT_VERSION}")
                endif()
            else()
                # Fallback to C++20 if version detection fails
                set(CMAKE_CXX_STANDARD 20)
                message(STATUS "Could not detect RDKit version, defaulting to C++20")
            endif()
        else()
            # Fallback to C++20 if version detection fails
            set(CMAKE_CXX_STANDARD 20)
            message(STATUS "Could not parse RDKit version, defaulting to C++20")
        endif()
    else()
        # Fallback to C++20 if no version header found
        set(CMAKE_CXX_STANDARD 20)
        message(STATUS "No RDKit version header found, defaulting to C++20")
    endif()
else()
    # Fallback to C++20 if RDKit not found initially
    set(CMAKE_CXX_STANDARD 20)
    message(STATUS "RDKit include directory not found initially, defaulting to C++20")
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)

# For static linking and portability
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(BUILD_SHARED_LIBS OFF)

# Set static linking flags for portable .so file
if(UNIX AND NOT APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc -static-libstdc++")
endif()

# Find required packages
find_package(PkgConfig REQUIRED)

# RDKit should already be found from above
if(RDKIT_INCLUDE_DIR AND RDKIT_GRAPHMOL_LIB)
    set(RDKit_FOUND TRUE)
    # Find all required RDKit libraries
    set(RDKIT_LIBRARIES)
    set(RDKIT_LIB_NAMES
        RDKitGraphMol
        RDKitMolTransforms
        RDKitForceFieldHelpers
        RDKitForceField
        RDKitDistGeom
        RDKitDistGeomHelpers
        RDKitEigenSolvers
        RDKitNumerics
        RDKitAlignment
        RDKitMMFFHelpers
        RDKitSmilesParse
        RDKitMolOps
        RDKitSubstructMatch
        RDKitQueryOps
        RDKitGeometry
        RDKitRDGeneral
    )
    
    foreach(lib_name ${RDKIT_LIB_NAMES})
        find_library(${lib_name}_LIB ${lib_name}
            PATHS ${RDKit_ROOT_DIR}/lib
                  /usr/lib
                  /usr/local/lib
            NO_DEFAULT_PATH)
        if(${lib_name}_LIB)
            list(APPEND RDKIT_LIBRARIES ${${lib_name}_LIB})
            message(STATUS "Found ${lib_name}: ${${lib_name}_LIB}")
        endif()
    endforeach()
else()
    # Fallback to package-based search
    find_package(RDKit QUIET)
endif()

if(NOT RDKit_FOUND)
    message(FATAL_ERROR "RDKit not found. Please ensure RDKit is installed and RDBASE or CONDA_PREFIX is set.")
endif()

# Find PyBind11
find_package(pybind11 QUIET)
if(NOT pybind11_FOUND)
    # Try to find pybind11 in conda environment
    find_path(PYBIND11_INCLUDE_DIR pybind11/pybind11.h
        PATHS ${RDKit_ROOT_DIR}/include
              /usr/include
              /usr/local/include
        NO_DEFAULT_PATH)
    
    if(PYBIND11_INCLUDE_DIR)
        set(pybind11_FOUND TRUE)
        message(STATUS "Found PyBind11 headers: ${PYBIND11_INCLUDE_DIR}")
        
        # Import pybind11 targets manually
        add_library(pybind11::headers INTERFACE IMPORTED)
        set_target_properties(pybind11::headers PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${PYBIND11_INCLUDE_DIR}")
        
        # Find Python
        find_package(Python COMPONENTS Interpreter Development REQUIRED)
        
        # Create pybind11::module target
        add_library(pybind11::module INTERFACE IMPORTED)
        set_target_properties(pybind11::module PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${PYBIND11_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "Python::Module")
        
        # Set compile definitions for pybind11 module
        target_compile_definitions(pybind11::module INTERFACE
            VERSION_INFO="dev")
    endif()
endif()

if(NOT pybind11_FOUND)
    message(FATAL_ERROR "PyBind11 not found. Please ensure PyBind11 is installed.")
endif()

# Find Boost (required by RDKit)
find_package(Boost REQUIRED COMPONENTS system filesystem serialization)

# Find Eigen3 (required by RDKit)
find_package(Eigen3 QUIET)
if(NOT Eigen3_FOUND)
    find_path(EIGEN3_INCLUDE_DIR Eigen/Core
        PATHS ${RDKit_ROOT_DIR}/include/eigen3
              ${RDKit_ROOT_DIR}/include
              /usr/include/eigen3
              /usr/local/include/eigen3
        NO_DEFAULT_PATH)
    
    if(EIGEN3_INCLUDE_DIR)
        set(Eigen3_FOUND TRUE)
        message(STATUS "Found Eigen3 headers: ${EIGEN3_INCLUDE_DIR}")
    endif()
endif()

# Create the library
add_library(stochastic_sampling_lib STATIC
    common_types.cpp
    accelerated_rmsd.cpp
    stochastic_sampling.cpp
)

# Set target properties for the library
target_include_directories(stochastic_sampling_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link RDKit to the library
if(RDKIT_INCLUDE_DIR)
    target_include_directories(stochastic_sampling_lib PUBLIC ${RDKIT_INCLUDE_DIR})
    target_link_libraries(stochastic_sampling_lib PUBLIC ${RDKIT_LIBRARIES})
else()
    target_link_libraries(stochastic_sampling_lib PUBLIC RDKit::GraphMol RDKit::MolTransforms RDKit::ForceFieldHelpers)
endif()

# Link Boost
target_link_libraries(stochastic_sampling_lib PUBLIC Boost::system Boost::filesystem Boost::serialization)

# Link Eigen3
if(EIGEN3_INCLUDE_DIR)
    target_include_directories(stochastic_sampling_lib PUBLIC ${EIGEN3_INCLUDE_DIR})
elseif(TARGET Eigen3::Eigen)
    target_link_libraries(stochastic_sampling_lib PUBLIC Eigen3::Eigen)
endif()

# Create the Python module
pybind11_add_module(stochastic_sampling_combined 
    pybind_interface.cpp
)

# Link the static library to the Python module
target_link_libraries(stochastic_sampling_combined PRIVATE stochastic_sampling_lib)

# Properties for the Python module
target_compile_definitions(stochastic_sampling_combined PRIVATE VERSION_INFO="1.0.0")

# Set properties for portable .so file
set_target_properties(stochastic_sampling_combined PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    VISIBILITY_INLINES_HIDDEN ON
)

# Set RPATH to find RDKit libraries
if(RDKIT_LIBRARIES)
    get_filename_component(RDKIT_LIB_DIR "${RDKIT_LIBRARIES}" DIRECTORY)
    set_target_properties(stochastic_sampling_combined PROPERTIES
        INSTALL_RPATH "${RDKIT_LIB_DIR}"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
endif()

# Force static linking of standard libraries
if(UNIX AND NOT APPLE)
    target_link_options(stochastic_sampling_combined PRIVATE
        -static-libgcc
        -static-libstdc++
    )
endif()


# Print configuration summary
message(STATUS "")
message(STATUS "=== Configuration Summary ===")
message(STATUS "RDKit found: ${RDKit_FOUND}")
if(RDKIT_INCLUDE_DIR)
    message(STATUS "RDKit include: ${RDKIT_INCLUDE_DIR}")
    message(STATUS "RDKit libraries: ${RDKIT_LIBRARIES}")
endif()
message(STATUS "PyBind11 found: ${pybind11_FOUND}")
if(PYBIND11_INCLUDE_DIR)
    message(STATUS "PyBind11 include: ${PYBIND11_INCLUDE_DIR}")
endif()
message(STATUS "Boost found: ${Boost_FOUND}")
message(STATUS "Eigen3 found: ${Eigen3_FOUND}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Static linking: ON (for portability)")
message(STATUS "=============================")
message(STATUS "")

# Add custom target for building and testing
add_custom_target(build_and_test
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target stochastic_sampling_combined
    COMMAND ${CMAKE_COMMAND} -E echo "Module built successfully"
    COMMAND ${CMAKE_COMMAND} -E echo "To test, run: python test_combined.py -s 'CCO' -d -nconfs 10"
    DEPENDS stochastic_sampling_combined
)

# Create portable deployment package
add_custom_target(create_portable_package
    COMMAND ${CMAKE_COMMAND} -E make_directory portable_deployment
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:stochastic_sampling_combined> portable_deployment/
    COMMAND ${CMAKE_COMMAND} -E echo "Copying RDKit libraries to portable_deployment..."
    COMMAND ${CMAKE_COMMAND} -E copy ${RDKIT_LIBRARIES} portable_deployment/
    COMMAND ${CMAKE_COMMAND} -E echo "Portable deployment created in: portable_deployment/"
    COMMAND ${CMAKE_COMMAND} -E echo "To deploy:"
    COMMAND ${CMAKE_COMMAND} -E echo "1. Copy the entire portable_deployment/ directory to target system"
    COMMAND ${CMAKE_COMMAND} -E echo "2. Set LD_LIBRARY_PATH=\$PWD/portable_deployment or add to system library path"
    COMMAND ${CMAKE_COMMAND} -E echo "3. Run your Python script"
    DEPENDS stochastic_sampling_combined
)