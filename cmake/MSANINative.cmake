include_guard(GLOBAL)

# Packaging builds use a private, static RDKit SDK. Keep the default shared
# mode for local development until the static SDK has been provisioned.
set(MSANI_RDKIT_LINKAGE "shared" CACHE STRING
    "RDKit linkage mode: shared (development) or static (redistributable wheel)")
set_property(CACHE MSANI_RDKIT_LINKAGE PROPERTY STRINGS shared static)

function(msani_validate_rdkit_linkage)
    if(NOT MSANI_RDKIT_LINKAGE STREQUAL "shared" AND
       NOT MSANI_RDKIT_LINKAGE STREQUAL "static")
        message(FATAL_ERROR
            "MSANI_RDKIT_LINKAGE must be 'shared' or 'static', got "
            "'${MSANI_RDKIT_LINKAGE}'")
    endif()

    if(NOT MSANI_RDKIT_LINKAGE STREQUAL "static")
        return()
    endif()

    if(MSVC)
        message(FATAL_ERROR
            "The private static RDKit wheel path is currently implemented "
            "for Linux packaging only. Use MSANI_RDKIT_LINKAGE=shared for "
            "Windows development builds.")
    endif()

    foreach(_msani_rdkit_target IN LISTS ARGN)
        if(NOT TARGET ${_msani_rdkit_target})
            message(FATAL_ERROR
                "Static packaging requires imported target "
                "${_msani_rdkit_target}. Configure CMAKE_PREFIX_PATH with "
                "the pinned static RDKit SDK.")
        endif()

        get_target_property(_msani_rdkit_location
            ${_msani_rdkit_target} IMPORTED_LOCATION_RELEASE)
        if(NOT _msani_rdkit_location)
            get_target_property(_msani_rdkit_location
                ${_msani_rdkit_target} IMPORTED_LOCATION)
        endif()
        if(NOT _msani_rdkit_location MATCHES "\\.a$")
            message(FATAL_ERROR
                "Static packaging requires ${_msani_rdkit_target} to resolve "
                "to a static archive, but it resolves to "
                "'${_msani_rdkit_location}'. Build/install RDKit 2025.09.5 "
                "with BUILD_SHARED_LIBS=OFF and point CMAKE_PREFIX_PATH to it.")
        endif()
    endforeach()
endfunction()

function(msani_configure_native_static_target target)
    set_target_properties(${target} PROPERTIES
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN YES
        POSITION_INDEPENDENT_CODE ON
    )
endfunction()

function(msani_configure_native_extension target)
    set_target_properties(${target} PROPERTIES
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN YES
        SKIP_BUILD_RPATH TRUE
        BUILD_WITH_INSTALL_RPATH FALSE
        INSTALL_RPATH ""
        INSTALL_RPATH_USE_LINK_PATH FALSE
    )

    # Local static archives include MSANI's copied RDKit-adjacent code. Do not
    # export those symbols into the process-wide dynamic-linker namespace.
    if(UNIX AND NOT APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_link_options(${target} PRIVATE "-Wl,--exclude-libs,ALL")
    endif()
endfunction()
