
# Find or build g2o.
# Requires CMake 3.16+ and a single-configuration generator.
include_guard(GLOBAL)

# Build configuration
if(CMAKE_CONFIGURATION_TYPES)
    message(FATAL_ERROR "g2o fallback requires Ninja or Unix Makefiles")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Project directories
get_filename_component(PROJECT_ROOT
    "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE
)

set(EXTERNAL_BUILD_ROOT
    "${PROJECT_ROOT}/build/externals"
    CACHE PATH "External build root"
)

set(EXTERNAL_INSTALL_ROOT
    "${PROJECT_ROOT}/install/externals"
    CACHE PATH "External install root"
)

set(G2O_INSTALL_DIR
    "${EXTERNAL_INSTALL_ROOT}/g2o/${CMAKE_BUILD_TYPE}"
)

# 1. Check local installation
find_package(g2o CONFIG QUIET
    PATHS "${G2O_INSTALL_DIR}"
    NO_DEFAULT_PATH
)

# 2. Check system installation
if(NOT g2o_FOUND)
    find_package(g2o CONFIG QUIET)
endif()

# 3. Build if not found
if(NOT g2o_FOUND)
    message(STATUS "g2o not found. Building ${CMAKE_BUILD_TYPE}...")

    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env
            "BUILD_TYPE=${CMAKE_BUILD_TYPE}"
            "EXTERNAL_BUILD_ROOT=${EXTERNAL_BUILD_ROOT}"
            "EXTERNAL_INSTALL_ROOT=${EXTERNAL_INSTALL_ROOT}"
            "MARCH_NATIVE=${CG2O_G2O_MARCH_NATIVE}"
            bash "${PROJECT_ROOT}/scripts/build_g2o.sh"
        RESULT_VARIABLE G2O_BUILD_RESULT
    )

    if(NOT "${G2O_BUILD_RESULT}" STREQUAL "0")
        message(FATAL_ERROR "Failed to build g2o: ${G2O_BUILD_RESULT}")
    endif()

    # Load newly installed g2o
    unset(g2o_DIR CACHE)

    find_package(g2o CONFIG REQUIRED
        PATHS "${G2O_INSTALL_DIR}"
        NO_DEFAULT_PATH
    )
endif()

# 4. Verify targets
if(NOT TARGET g2o::core OR NOT TARGET g2o::stuff)
    message(FATAL_ERROR "g2o::core or g2o::stuff is missing")
endif()

set(G2O_LIBRARIES g2o::core g2o::stuff)

message(STATUS "g2o ready (${CMAKE_BUILD_TYPE})")
