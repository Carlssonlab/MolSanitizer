cmake_minimum_required(VERSION 3.16)
project(msani_confgen_cpp)

# Always find the Python interpreter from the build environment
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Find RDKit first to determine version for C++ standard
set(RDKit_ROOT_DIR "$ENV{RDBASE}")
if(NOT RDKit_ROOT_DIR)
    set(RDKit_ROOT_DIR "$ENV{CONDA_PREFIX}")
endif()

# Try to get Python's site-packages path
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import sys; import sysconfig; print(sysconfig.get_path('purelib'))"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
message(STATUS "Python site-packages: ${PYTHON_SITE_PACKAGES}")
message(STATUS "RDKit root directory: ${RDKit_ROOT_DIR}")
message(STATUS "Python executable: ${Python3_EXECUTABLE}")

# Try to detect RDKit version using Python to determine C++ standard
# Add explicit PYTHONPATH from conda environment to ensure rdkit can be found
set(ENV{PYTHONPATH} "${RDKit_ROOT_DIR}/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages:${PYTHON_SITE_PACKAGES}:$ENV{PYTHONPATH}")
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import sys; sys.path.insert(0, '${RDKit_ROOT_DIR}/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages'); sys.path.insert(0, '${PYTHON_SITE_PACKAGES}'); import rdkit; print(rdkit.__version__)"
    OUTPUT_VARIABLE RDKIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE RESULT
    ERROR_QUIET
)
if(RESULT EQUAL 0)
    message(STATUS "Detected RDKit version via Python: ${RDKIT_VERSION}")
    # Set C++ standard based on RDKit version
    # RDKit <= 2025.3.1: C++17, > 2025.3.1: C++20
    if(RDKIT_VERSION VERSION_LESS_EQUAL "2025.3.1")
        set(CMAKE_CXX_STANDARD 17)
        message(STATUS "Using C++17 for RDKit version ${RDKIT_VERSION}")
    else()
        set(CMAKE_CXX_STANDARD 20)
        message(STATUS "Using C++20 for RDKit version ${RDKIT_VERSION}")
    endif()
else()
    # Fallback to C++17 if Python detection fails
    set(CMAKE_CXX_STANDARD 17)
    message(STATUS "Could not detect RDKit version via Python, defaulting to C++17")
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Detect whether RDKit was built with thread-safe substructure search support
set(RDK_THREADSAFE_SSS OFF)
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import sys; sys.path.insert(0, '${RDKit_ROOT_DIR}/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages'); sys.path.insert(0, '${PYTHON_SITE_PACKAGES}'); from rdkit import rdBase; print('ON' if getattr(rdBase, '_multithreadedEnabled', False) else 'OFF')"
    OUTPUT_VARIABLE RDK_THREADSAFE_FLAG
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE RDK_THREADSAFE_RESULT
    ERROR_QUIET
)
if(RDK_THREADSAFE_RESULT EQUAL 0 AND RDK_THREADSAFE_FLAG STREQUAL "ON")
    set(RDK_THREADSAFE_SSS ON)
    message(STATUS "RDKit reports multithreaded embedding support: ON")
else()
    message(STATUS "RDKit reports multithreaded embedding support: OFF")
endif()

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

# Find RDKit - try multiple approaches
set(RDKit_ROOT_DIR "$ENV{RDBASE}")
if(NOT RDKit_ROOT_DIR)
    set(RDKit_ROOT_DIR "$ENV{CONDA_PREFIX}")
endif()

# Find RDKit - try multiple approaches
set(RDKit_ROOT_DIR "$ENV{RDBASE}")
if(NOT RDKit_ROOT_DIR)
    set(RDKit_ROOT_DIR "$ENV{CONDA_PREFIX}")
endif()

# Try to find RDKit directly without imported targets first
find_path(RDKIT_INCLUDE_DIR GraphMol/ROMol.h
    PATHS ${RDKit_ROOT_DIR}/include/rdkit
          ${RDKit_ROOT_DIR}/include
          /usr/include/rdkit
          /usr/local/include/rdkit
    NO_DEFAULT_PATH)

find_library(RDKIT_GRAPHMOL_LIB RDKitGraphMol
    PATHS ${RDKit_ROOT_DIR}/lib
          /usr/lib
          /usr/local/lib
    NO_DEFAULT_PATH)

if(RDKIT_INCLUDE_DIR AND RDKIT_GRAPHMOL_LIB)
    set(RDKit_FOUND TRUE)
    message(STATUS "Found RDKit headers: ${RDKIT_INCLUDE_DIR}")
    message(STATUS "Found RDKit GraphMol lib: ${RDKIT_GRAPHMOL_LIB}")
    
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
    initial_embedder.cpp
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
pybind11_add_module(msani_confgen_cpp 
    pybind_interface.cpp
)

# Link the static library to the Python module
target_link_libraries(msani_confgen_cpp PRIVATE stochastic_sampling_lib)

# Properties for the Python module
target_compile_definitions(msani_confgen_cpp PRIVATE VERSION_INFO="1.0.0")

if(RDK_THREADSAFE_SSS)
    target_compile_definitions(stochastic_sampling_lib PUBLIC RDK_BUILD_THREADSAFE_SSS)
    target_compile_definitions(msani_confgen_cpp PRIVATE RDK_BUILD_THREADSAFE_SSS)
endif()

# Set properties for portable .so file
set_target_properties(msani_confgen_cpp PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    VISIBILITY_INLINES_HIDDEN ON
)

# Set RPATH to find RDKit libraries
if(RDKIT_LIBRARIES)
    get_filename_component(RDKIT_LIB_DIR "${RDKIT_LIBRARIES}" DIRECTORY)
    set_target_properties(msani_confgen_cpp PROPERTIES
        INSTALL_RPATH "${RDKIT_LIB_DIR}"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
endif()

# Force static linking of standard libraries
if(UNIX AND NOT APPLE)
    target_link_options(msani_confgen_cpp PRIVATE
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
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target msani_confgen_cpp
    COMMAND ${CMAKE_COMMAND} -E echo "Module built successfully"
    COMMAND ${CMAKE_COMMAND} -E echo "To test, run: python test_combined.py -s 'CCO' -d -nconfs 10"
    DEPENDS msani_confgen_cpp
)

# Create portable deployment package
add_custom_target(create_portable_package
    COMMAND ${CMAKE_COMMAND} -E make_directory portable_deployment
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:msani_confgen_cpp> portable_deployment/
    COMMAND ${CMAKE_COMMAND} -E echo "Copying RDKit libraries to portable_deployment..."
    COMMAND ${CMAKE_COMMAND} -E copy ${RDKIT_LIBRARIES} portable_deployment/
    COMMAND ${CMAKE_COMMAND} -E echo "Portable deployment created in: portable_deployment/"
    COMMAND ${CMAKE_COMMAND} -E echo "To deploy:"
    COMMAND ${CMAKE_COMMAND} -E echo "1. Copy the entire portable_deployment/ directory to target system"
    COMMAND ${CMAKE_COMMAND} -E echo "2. Set LD_LIBRARY_PATH=\$PWD/portable_deployment or add to system library path"
    COMMAND ${CMAKE_COMMAND} -E echo "3. Run your Python script"
    DEPENDS msani_confgen_cpp
)

# Install the Python module
install(TARGETS msani_confgen_cpp
    LIBRARY DESTINATION .
    COMPONENT python
)
