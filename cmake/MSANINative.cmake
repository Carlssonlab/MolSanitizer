include_guard(GLOBAL)

# Packaging builds use a private, static RDKit SDK. Keep the default shared
# mode for local development until the static SDK has been provisioned.
set(MSANI_RDKIT_LINKAGE "shared" CACHE STRING
    "RDKit linkage mode: shared (development) or static (redistributable wheel)")
set_property(CACHE MSANI_RDKIT_LINKAGE PROPERTY STRINGS shared static)

function(msani_select_rdkit_targets output_variable)
    if(NOT MSANI_RDKIT_LINKAGE STREQUAL "shared" AND
       NOT MSANI_RDKIT_LINKAGE STREQUAL "static")
        message(FATAL_ERROR
            "MSANI_RDKIT_LINKAGE must be 'shared' or 'static', got "
            "'${MSANI_RDKIT_LINKAGE}'")
    endif()

    set(_msani_selected_targets)
    foreach(_msani_rdkit_component IN LISTS ARGN)
        if(MSANI_RDKIT_LINKAGE STREQUAL "static")
            set(_msani_rdkit_target "RDKit::${_msani_rdkit_component}_static")
        else()
            set(_msani_rdkit_target "RDKit::${_msani_rdkit_component}")
        endif()

        if(NOT TARGET ${_msani_rdkit_target})
            message(FATAL_ERROR
                "RDKit linkage mode '${MSANI_RDKIT_LINKAGE}' requires target "
                "${_msani_rdkit_target}. Point CMAKE_PREFIX_PATH or RDBASE "
                "at the pinned RDKit SDK.")
        endif()

        if(MSANI_RDKIT_LINKAGE STREQUAL "static")
            get_target_property(_msani_rdkit_type ${_msani_rdkit_target} TYPE)
            if(NOT _msani_rdkit_type STREQUAL "STATIC_LIBRARY")
                message(FATAL_ERROR
                    "${_msani_rdkit_target} is ${_msani_rdkit_type}, not a "
                    "static RDKit library. Rebuild the selected RDKit SDK with "
                    "RDK_BUILD_STATIC_LIBS_ONLY=ON.")
            endif()
        endif()

        list(APPEND _msani_selected_targets ${_msani_rdkit_target})
    endforeach()

    set(${output_variable} "${_msani_selected_targets}" PARENT_SCOPE)
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
