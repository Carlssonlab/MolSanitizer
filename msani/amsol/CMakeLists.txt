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

install(TARGETS _amsolcpp
    LIBRARY DESTINATION amsolcpp COMPONENT python
    RUNTIME DESTINATION amsolcpp COMPONENT python
)
install(FILES python/amsolcpp/__init__.py
    DESTINATION amsolcpp
    COMPONENT python
)
