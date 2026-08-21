cmake_minimum_required(VERSION 3.16)
project(msani_confgen_cpp)

# Always find the Python interpreter from the build environment
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Find RDKit root (Conda on Windows uses Library/*)
set(RDKit_ROOT_DIR "$ENV{RDBASE}")
if(NOT RDKit_ROOT_DIR)
    set(RDKit_ROOT_DIR "$ENV{CONDA_PREFIX}")
endif()

# Fallback: derive environment root from Python's sys.prefix if still empty
if(NOT RDKit_ROOT_DIR)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import sys; print(sys.prefix)"
        OUTPUT_VARIABLE _PY_SYS_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(_PY_SYS_PREFIX)
        set(RDKit_ROOT_DIR "${_PY_SYS_PREFIX}")
    endif()
endif()

# Help CMake find Conda-provided CMake packages (RDKitConfig.cmake, pybind11Config.cmake)
if(RDKit_ROOT_DIR)
    list(PREPEND CMAKE_PREFIX_PATH
        "${RDKit_ROOT_DIR}"
        "${RDKit_ROOT_DIR}/Library"
        "${RDKit_ROOT_DIR}/Library/lib/cmake"
        "${RDKit_ROOT_DIR}/lib/cmake"
    )
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
if(WIN32)
    set(PYTHON_PATH_SEP ";")
    set(RDKIT_SITE_PACKAGES "${RDKit_ROOT_DIR}/Lib/site-packages")
    # Ensure RDKit DLLs are in PATH for the import to work
    set(ENV{PATH} "${RDKit_ROOT_DIR}/Library/bin;$ENV{PATH}")
else()
    set(PYTHON_PATH_SEP ":")
    set(RDKIT_SITE_PACKAGES "${RDKit_ROOT_DIR}/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages")
endif()

set(ENV{PYTHONPATH} "${RDKIT_SITE_PACKAGES}${PYTHON_PATH_SEP}${PYTHON_SITE_PACKAGES}${PYTHON_PATH_SEP}$ENV{PYTHONPATH}")
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import sys; sys.path.insert(0, '${RDKIT_SITE_PACKAGES}'); sys.path.insert(0, '${PYTHON_SITE_PACKAGES}'); import rdkit; print(rdkit.__version__)"
    OUTPUT_VARIABLE RDKIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE RESULT
    ERROR_QUIET
)
if(RESULT EQUAL 0)
    message(STATUS "Detected RDKit version via Python: ${RDKIT_VERSION}")
    # Set C++ standard based on RDKit version
    # RDKit < 2025.3.1: C++17, >= 2025.3.1: C++20
    if(RDKIT_VERSION VERSION_LESS "2025.3.1")
        set(CMAKE_CXX_STANDARD 17)
        message(STATUS "Using C++17 for RDKit version ${RDKIT_VERSION}")
    else()
        set(CMAKE_CXX_STANDARD 20)
        message(STATUS "Using C++20 for RDKit version ${RDKIT_VERSION}")
    endif()
else()
    # Fallback to C++20 if Python detection fails
    set(CMAKE_CXX_STANDARD 20)
    message(STATUS "Could not detect RDKit version via Python, defaulting to C++20")
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Detect whether RDKit was built with thread-safe substructure search support
set(RDK_THREADSAFE_SSS OFF)
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import sys; sys.path.insert(0, '${RDKIT_SITE_PACKAGES}'); sys.path.insert(0, '${PYTHON_SITE_PACKAGES}'); from rdkit import rdBase; print('ON' if getattr(rdBase, '_multithreadedEnabled', False) else 'OFF')"
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

# All local static libraries must be usable from Python extension modules.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Platform-specific compiler flags for optimization
if(MSVC)
    # Microsoft Visual C++ - Remove all optimization flags that might interfere
    # Clear any inherited flags
    string(REPLACE "/O2" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
    string(REPLACE "/Ob2" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /Ob2 /DNDEBUG" CACHE STRING "" FORCE)
    # Use target-specific flags instead of global add_compile_options to avoid conflicts
    set(MSVC_COMPILE_FLAGS /W3 /EHsc /bigobj /Zc:__cplusplus)
    # Define BOOST_NO_EXCEPTIONS=0 to ensure Boost exception handling is enabled
    add_compile_definitions(BOOST_NO_EXCEPTIONS=0)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # GCC or Clang
    add_compile_options(-O3 -Wall)
    set(MSVC_COMPILE_FLAGS "")
endif()

# Some build frontends inject '-O3' unconditionally; remove it for MSVC to avoid cl warnings
if(MSVC)
    foreach(flag_var
        CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE CMAKE_C_FLAGS_RELWITHDEBINFO CMAKE_C_FLAGS_MINSIZEREL
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_RELWITHDEBINFO CMAKE_CXX_FLAGS_MINSIZEREL)
        if(DEFINED ${flag_var})
            string(REPLACE "-O3" "" ${flag_var} "${${flag_var}}")
            set(${flag_var} "${${flag_var}}" CACHE STRING "" FORCE)
        endif()
    endforeach()
endif()

# Find required packages
find_package(PkgConfig QUIET)

## RDKit discovery
find_package(RDKit CONFIG QUIET)

# Try to find RDKit directly without imported targets first
if(NOT RDKit_FOUND)
    find_path(RDKIT_INCLUDE_DIR GraphMol/ROMol.h
        PATHS
            ${RDKit_ROOT_DIR}/include/rdkit
            ${RDKit_ROOT_DIR}/include
            ${RDKit_ROOT_DIR}/Library/include/rdkit
            ${RDKit_ROOT_DIR}/Library/include
            /usr/include/rdkit
            /usr/local/include/rdkit
    )

    find_library(RDKIT_GRAPHMOL_LIB NAMES RDKitGraphMol
        PATHS
            ${RDKit_ROOT_DIR}/lib
            ${RDKit_ROOT_DIR}/Library/lib
            /usr/lib
            /usr/local/lib
    )
endif()

if(RDKit_FOUND OR (RDKIT_INCLUDE_DIR AND RDKIT_GRAPHMOL_LIB))
    set(RDKit_FOUND TRUE)
    if(RDKIT_INCLUDE_DIR)
        message(STATUS "Found RDKit headers: ${RDKIT_INCLUDE_DIR}")
    endif()
    if(RDKIT_GRAPHMOL_LIB)
        message(STATUS "Found RDKit GraphMol lib: ${RDKIT_GRAPHMOL_LIB}")
    endif()
    
    # Find all required RDKit libraries
    set(RDKIT_LIBRARIES)
    set(RDKIT_LIB_NAMES
        RDKitGraphMol
        RDKitMolTransforms
        RDKitForceFieldHelpers
        RDKitForceField
        RDKitDistGeomHelpers
        RDKitDistGeometry
        RDKitEigenSolvers
        RDKitAlignment
        RDKitSmilesParse
        RDKitSubstructMatch
        RDKitRDGeneral
        RDKitRDGeometryLib
    )
    
    if(NOT RDKit_FOUND)
        foreach(lib_name ${RDKIT_LIB_NAMES})
            find_library(${lib_name}_LIB NAMES ${lib_name}
                PATHS
                    ${RDKit_ROOT_DIR}/lib
                    ${RDKit_ROOT_DIR}/Library/lib
                    /usr/lib
                    /usr/local/lib
            )
            if(${lib_name}_LIB)
                list(APPEND RDKIT_LIBRARIES ${${lib_name}_LIB})
                message(STATUS "Found ${lib_name}: ${${lib_name}_LIB}")
            endif()
        endforeach()
    endif()
endif()

if(NOT RDKit_FOUND)
    message(FATAL_ERROR "RDKit not found. Please ensure RDKit is installed and RDBASE or CONDA_PREFIX is set.")
endif()

# Wheel builds use a pinned static RDKit SDK; the shared/manual lookup path is
# retained for local development only.
msani_select_rdkit_targets(MSANI_CONFGEN_RDKIT_TARGETS
    RDGeneral RDGeometryLib SubstructMatch SmilesParse GraphMol EigenSolvers
    Alignment ForceField ForceFieldHelpers DistGeometry DistGeomHelpers
    MolTransforms
)

# Find PyBind11
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    # Try to find pybind11 in conda environment
    find_path(PYBIND11_INCLUDE_DIR pybind11/pybind11.h
        PATHS
            ${RDKit_ROOT_DIR}/include
            ${RDKit_ROOT_DIR}/Library/include
            /usr/include
            /usr/local/include
    )
    
    if(PYBIND11_INCLUDE_DIR)
        set(pybind11_FOUND TRUE)
        message(STATUS "Found PyBind11 headers: ${PYBIND11_INCLUDE_DIR}")
        
        # Import pybind11 targets manually
        add_library(pybind11::headers INTERFACE IMPORTED)
        set_target_properties(pybind11::headers PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${PYBIND11_INCLUDE_DIR}")
        
    # Find Python (already found as Python3)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
        
        # Create pybind11::module target
        add_library(pybind11::module INTERFACE IMPORTED)
        set_target_properties(pybind11::module PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${PYBIND11_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "Python3::Module")
        
        # Set compile definitions for pybind11 module
        target_compile_definitions(pybind11::module INTERFACE
            VERSION_INFO="dev")
    endif()
endif()

if(NOT pybind11_FOUND)
    message(FATAL_ERROR "PyBind11 not found. Please ensure PyBind11 is installed.")
endif()

# Find Boost (required by RDKit)
find_package(Boost REQUIRED)

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
    $<$<PLATFORM_ID:Windows>:boost_throw_stub.cpp>
    $<$<PLATFORM_ID:Windows>:rdkit_romol_destroy_stub.cpp>
)
msani_configure_native_static_target(stochastic_sampling_lib)

# Apply MSVC-specific compile flags to this target
if(MSVC)
    target_compile_options(stochastic_sampling_lib PRIVATE ${MSVC_COMPILE_FLAGS})
    target_compile_definitions(stochastic_sampling_lib PRIVATE RDKIT_DYN_LINK)
endif()

# Set target properties for the library
target_include_directories(stochastic_sampling_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link RDKit to the library
if(RDKit_FOUND)
    # Use imported RDKit targets (provides both includes and libs automatically)
    # Link in dependency order: low-level libraries first
    target_link_libraries(stochastic_sampling_lib PUBLIC ${MSANI_CONFGEN_RDKIT_TARGETS})
    # On Windows with MSVC, RDKit is built as DLLs
    # Add the DLL directory to runtime path
    if(MSVC AND RDKit_ROOT_DIR)
        set_target_properties(stochastic_sampling_lib PROPERTIES
            VS_DEBUGGER_ENVIRONMENT "PATH=${RDKit_ROOT_DIR}/Library/bin;$ENV{PATH}"
        )
    endif()
else()
    # Fallback: manual include and lib linking
    if(RDKIT_INCLUDE_DIR)
        target_include_directories(stochastic_sampling_lib PUBLIC ${RDKIT_INCLUDE_DIR})
    endif()
    target_link_libraries(stochastic_sampling_lib PUBLIC ${RDKIT_LIBRARIES})
endif()

# We don't link Boost directly since we only use header-only parts.
# RDKit DLLs handle their own transitive Boost dependencies.

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

# Apply MSVC-specific compile flags to this target
if(MSVC)
    target_compile_options(msani_confgen_cpp PRIVATE ${MSVC_COMPILE_FLAGS})
    target_compile_definitions(msani_confgen_cpp PRIVATE RDKIT_DYN_LINK)
    # Enable verbose library search output from the linker
    target_link_options(msani_confgen_cpp PRIVATE /VERBOSE:LIB)
endif()

# Link the static library to the Python module
target_link_libraries(msani_confgen_cpp PRIVATE stochastic_sampling_lib)

# On Windows/MSVC, explicitly link RDKit and Boost to the extension as well
if(RDKit_FOUND)
    target_link_libraries(msani_confgen_cpp PRIVATE ${MSANI_CONFGEN_RDKIT_TARGETS})
else()
    if(RDKIT_LIBRARIES)
        target_link_libraries(msani_confgen_cpp PRIVATE ${RDKIT_LIBRARIES})
    endif()
endif()

# Transitive linkage for Boost is handled by RDKit if needed.

# Properties for the Python module
target_compile_definitions(msani_confgen_cpp PRIVATE VERSION_INFO="1.0.0")

if(RDK_THREADSAFE_SSS)
    target_compile_definitions(stochastic_sampling_lib PUBLIC RDK_BUILD_THREADSAFE_SSS)
    target_compile_definitions(msani_confgen_cpp PRIVATE RDK_BUILD_THREADSAFE_SSS)
endif()

msani_configure_native_extension(msani_confgen_cpp)

# On Windows with MSVC, ensure RDKit DLLs are in the PATH at runtime
# The conda environment already has them in Library/bin
if(MSVC)
    # Copy RDKit DLLs next to the built extension for reliable loading
    if(DEFINED RDKit_ROOT_DIR)
        file(GLOB RDKIT_DLLS "${RDKit_ROOT_DIR}/Library/bin/RDKit*.dll")
        foreach(dll ${RDKIT_DLLS})
            add_custom_command(TARGET msani_confgen_cpp POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${dll}"
                    $<TARGET_FILE_DIR:msani_confgen_cpp>
            )
        endforeach()
    endif()
    # Add a note in the configuration output
    message(STATUS "")
    message(STATUS "Windows DLL Note: RDKit DLLs will be loaded from:")
    message(STATUS "  ${RDKit_ROOT_DIR}/Library/bin")
    message(STATUS "Ensure this directory is in your PATH when using the module.")
    message(STATUS "")
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
message(STATUS "RDKit linkage mode: ${MSANI_RDKIT_LINKAGE}")
message(STATUS "=============================")
message(STATUS "")

# Print link libraries for diagnostics
get_target_property(_msani_confgen_libs msani_confgen_cpp LINK_LIBRARIES)
message(STATUS "msani_confgen_cpp link libs: ${_msani_confgen_libs}")

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

# Install the .pyi stub file alongside the compiled module
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/../cpp/msani_confgen_cpp.pyi
    DESTINATION .
    COMPONENT python
)
