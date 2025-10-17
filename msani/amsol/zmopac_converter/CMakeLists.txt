cmake_minimum_required(VERSION 3.15)
project(ZmatMOPACWriter VERSION 1.0.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Platform detection
if(WIN32)
    set(PLATFORM_NAME "windows")
elseif(APPLE)
    set(PLATFORM_NAME "macos")
elseif(UNIX)
    set(PLATFORM_NAME "linux")
endif()

# Add executable
add_executable(zmopac_converter ZmatMOPACWriter.cpp)

# Link libraries
if(UNIX AND NOT APPLE)
    # Linux
    target_link_libraries(zmopac_converter m)
    # Static linking for better portability
    set_target_properties(zmopac_converter PROPERTIES
        LINK_FLAGS "-static-libgcc -static-libstdc++"
    )
elseif(APPLE)
    # macOS - link math library
    target_link_libraries(zmopac_converter m)
    # Universal binary support
    if(CMAKE_OSX_ARCHITECTURES)
        set_target_properties(zmopac_converter PROPERTIES
            OSX_ARCHITECTURES "${CMAKE_OSX_ARCHITECTURES}"
        )
    endif()
elseif(WIN32)
    # Windows - static runtime for portability
    if(MSVC)
        set_target_properties(zmopac_converter PROPERTIES
            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>"
        )
    endif()
endif()

# Compiler-specific options
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(zmopac_converter PRIVATE 
        -Wall -Wextra -O3 -DNDEBUG
        $<$<PLATFORM_ID:Linux>:-march=x86-64>
    )
elseif(MSVC)
    target_compile_options(zmopac_converter PRIVATE 
        /W4 /O2 /DNDEBUG /MT
    )
endif()

# Set output directory with platform-specific naming
set_target_properties(zmopac_converter PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Platform-specific executable naming
if(WIN32)
    set_target_properties(zmopac_converter PROPERTIES
        OUTPUT_NAME "zmopac_converter"
        SUFFIX ".exe"
    )
else()
    set_target_properties(zmopac_converter PROPERTIES
        OUTPUT_NAME "zmopac_converter"
    )
endif()

# Install target
install(TARGETS zmopac_converter
    RUNTIME DESTINATION bin
)
