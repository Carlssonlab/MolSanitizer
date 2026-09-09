cmake_minimum_required(VERSION 3.16)
project(msani_stereoisomers)

# Always find the Python interpreter from the build environment
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)

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
    message(WARNING "Could not detect RDKit version via Python, defaulting to C++20")
    set(CMAKE_CXX_STANDARD 20)
endif()

# Local static libraries must be usable from Python extension modules.
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
    # Leave BOOST_NO_EXCEPTIONS undefined: Boost checks its presence, not value.
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

#############################
# Dependency discovery
#############################

# Detect RDKit root from environment (Conda installs on Windows use Library/*)
set(RDKit_ROOT_DIR "$ENV{RDBASE}")
if(NOT RDKit_ROOT_DIR)
    set(RDKit_ROOT_DIR "$ENV{CONDA_PREFIX}")
endif()
message(STATUS "RDKit root directory: ${RDKit_ROOT_DIR}")

# Fallback: derive root from Python sys.prefix if env vars are missing
if(NOT RDKit_ROOT_DIR)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import sys; print(sys.prefix)"
        OUTPUT_VARIABLE _PY_SYS_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(_PY_SYS_PREFIX)
        set(RDKit_ROOT_DIR "${_PY_SYS_PREFIX}")
        message(STATUS "RDKit root derived from Python sys.prefix: ${RDKit_ROOT_DIR}")
    endif()
endif()

# Help CMake find Conda-provided CMake packages (RDKit, pybind11, etc.)
if(RDKit_ROOT_DIR)
    list(PREPEND CMAKE_PREFIX_PATH
        "${RDKit_ROOT_DIR}"
        "${RDKit_ROOT_DIR}/Library"
        "${RDKit_ROOT_DIR}/Library/lib/cmake"
        "${RDKit_ROOT_DIR}/lib/cmake"
    )
endif()

# Find Python and pybind11 (pybind11 provides a CMake config in Conda)
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# Prefer imported RDKit targets via CMake config (Conda installs ship RDKitConfig.cmake)
find_package(RDKit CONFIG QUIET)

# Fallback: manual header/lib discovery
if(NOT RDKit_FOUND)
    # Find RDKit include directory
    find_path(RDKIT_INCLUDE_DIR GraphMol/ROMol.h
        PATHS
            ${RDKit_ROOT_DIR}/include/rdkit
            ${RDKit_ROOT_DIR}/include
            ${RDKit_ROOT_DIR}/Library/include/rdkit
            ${RDKit_ROOT_DIR}/Library/include
            /usr/include/rdkit
            /usr/local/include/rdkit
    )

    # Find RDKit libraries (Windows .lib in Library/lib; Unix in lib)
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
        else()
            message(WARNING "Could not find ${lib_name}")
        endif()
    endforeach()
endif()

msani_select_rdkit_targets(MSANI_STEREO_RDKIT_TARGETS
    RDGeneral SmilesParse GraphMol ForceField ForceFieldHelpers DistGeomHelpers
)

# Find Boost (required by RDKit)
find_package(Boost REQUIRED)

# Create static library for stereoisomer enumeration
add_library(stereoisomer_lib STATIC
    EnumerateStereoisomers.cpp
    Flippers.cpp
)
msani_configure_native_static_target(stereoisomer_lib)
if(TARGET Boost::headers)
    target_link_libraries(stereoisomer_lib PUBLIC Boost::headers)
elseif(TARGET Boost::boost)
    target_link_libraries(stereoisomer_lib PUBLIC Boost::boost)
else()
    target_include_directories(stereoisomer_lib PUBLIC ${Boost_INCLUDE_DIRS})
endif()

# Apply MSVC-specific compile flags to this target
if(MSVC)
    target_compile_options(stereoisomer_lib PRIVATE ${MSVC_COMPILE_FLAGS})
    target_compile_definitions(stereoisomer_lib PRIVATE
        $<$<STREQUAL:${MSANI_RDKIT_LINKAGE},shared>:RDKIT_DYN_LINK>)
endif()

target_include_directories(stereoisomer_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link RDKit
if(RDKit_FOUND)
    # Link in dependency order: low-level libraries first
    target_link_libraries(stereoisomer_lib PUBLIC ${MSANI_STEREO_RDKIT_TARGETS})
    # On Windows with MSVC, RDKit is built as DLLs
    # Add the DLL directory to runtime path
    if(MSVC AND RDKit_ROOT_DIR)
        set_target_properties(stereoisomer_lib PROPERTIES
            VS_DEBUGGER_ENVIRONMENT "PATH=${RDKit_ROOT_DIR}/Library/bin;$ENV{PATH}"
        )
    endif()
else()
    if(RDKIT_INCLUDE_DIR)
        target_include_directories(stereoisomer_lib PUBLIC ${RDKIT_INCLUDE_DIR})
    endif()
    target_link_libraries(stereoisomer_lib PUBLIC ${RDKIT_LIBRARIES})
endif()

# We only use header-only Boost parts; RDKit DLLs handle transitive Boost deps.

# Create the Python module
pybind11_add_module(msani_stereoisomers
    ../bindings/stereoisomers.cpp
)

# Apply MSVC-specific compile flags to this target
if(MSVC)
    target_compile_options(msani_stereoisomers PRIVATE ${MSVC_COMPILE_FLAGS})
    target_compile_definitions(msani_stereoisomers PRIVATE
        $<$<STREQUAL:${MSANI_RDKIT_LINKAGE},shared>:RDKIT_DYN_LINK>)
    target_link_options(msani_stereoisomers PRIVATE /VERBOSE:LIB)
endif()

target_link_libraries(msani_stereoisomers PRIVATE
    stereoisomer_lib
)

# Explicitly link RDKit and Boost to the extension on MSVC (ensure import libs are included)
if(RDKit_FOUND)
    target_link_libraries(msani_stereoisomers PRIVATE ${MSANI_STEREO_RDKIT_TARGETS})
else()
    if(RDKIT_LIBRARIES)
        target_link_libraries(msani_stereoisomers PRIVATE ${RDKIT_LIBRARIES})
    endif()
endif()

# Transitive Boost linkage handled by RDKit.

# On Windows, copy RDKit DLLs next to the built extension for reliable loading
if(MSVC)
    if(DEFINED RDKit_ROOT_DIR)
        file(GLOB RDKIT_DLLS "${RDKit_ROOT_DIR}/Library/bin/RDKit*.dll")
        foreach(dll ${RDKIT_DLLS})
            add_custom_command(TARGET msani_stereoisomers POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${dll}"
                    $<TARGET_FILE_DIR:msani_stereoisomers>
            )
        endforeach()
    endif()
endif()

target_compile_definitions(msani_stereoisomers PRIVATE VERSION_INFO="0.4.0")
msani_configure_native_extension(msani_stereoisomers)

# On Windows with MSVC, ensure RDKit DLLs are in the PATH at runtime
# The conda environment already has them in Library/bin
if(MSVC)
    # Add a note in the configuration output
    message(STATUS "")
    message(STATUS "Windows DLL Note: RDKit DLLs will be loaded from:")
    message(STATUS "  ${RDKit_ROOT_DIR}/Library/bin")
    message(STATUS "Ensure this directory is in your PATH when using the module.")
    message(STATUS "")
endif()

# Install the Python module
install(TARGETS msani_stereoisomers
    LIBRARY DESTINATION .
    COMPONENT python
)

# Install the .pyi stub file alongside the compiled module
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/msani_stereoisomers.pyi
    DESTINATION .
    COMPONENT python
)

# Print link libraries for diagnostics
get_target_property(_msani_stereo_libs msani_stereoisomers LINK_LIBRARIES)
message(STATUS "msani_stereoisomers link libs: ${_msani_stereo_libs}")
