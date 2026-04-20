# FindCG2O.cmake
#
# Usage:
#   find_package(CG2O REQUIRED COMPONENTS core)
#   find_package(CG2O REQUIRED COMPONENTS core umfpack)
#   find_package(CG2O REQUIRED COMPONENTS core pardiso)
#   find_package(CG2O REQUIRED COMPONENTS cg2o)   # optional umbrella
#
# Result variables:
#   CG2O_FOUND
#   CG2O_INCLUDE_DIR
#   CG2O_LIBRARIES
#   CG2O_VERSION
#
# Imported targets:
#   CG2O::core
#   CG2O::umfpack
#   CG2O::pardiso
#   CG2O::cg2o
#
# Optional hint:
#   set(CG2O_ROOT "/path/to/install/prefix")

include(FindPackageHandleStandardArgs)

# ------------------------------------------------------------
# Hints
# ------------------------------------------------------------

set(_CG2O_HINTS)
if(CG2O_ROOT)
  list(APPEND _CG2O_HINTS "${CG2O_ROOT}")
endif()
if(DEFINED ENV{CG2O_ROOT})
  list(APPEND _CG2O_HINTS "$ENV{CG2O_ROOT}")
endif()

# ------------------------------------------------------------
# Includes
# ------------------------------------------------------------

find_path(CG2O_INCLUDE_DIR
  NAMES cg2o/core/base_fixed_sized_edge_eq.h
  HINTS ${_CG2O_HINTS}
  PATH_SUFFIXES include
)

# ------------------------------------------------------------
# Core library
# ------------------------------------------------------------

find_library(CG2O_CORE_LIBRARY
  NAMES
    cg2o_core
    core
  HINTS ${_CG2O_HINTS}
  PATH_SUFFIXES lib lib64
)

# ------------------------------------------------------------
# Optional umbrella library
# ------------------------------------------------------------

find_library(CG2O_CG2O_LIBRARY
  NAMES
    cg2o
    cg2o_cg2o
  HINTS ${_CG2O_HINTS}
  PATH_SUFFIXES lib lib64
)

# ------------------------------------------------------------
# Optional UMFPACK solver library
# ------------------------------------------------------------

find_library(CG2O_UMFPACK_LIBRARY
  NAMES
    solver_eigen_umfpack
    cg2o_solver_eigen_umfpack
    cg2o_solver_umfpack
  HINTS ${_CG2O_HINTS}
  PATH_SUFFIXES lib lib64
)

# ------------------------------------------------------------
# Optional Pardiso solver library
# ------------------------------------------------------------

find_library(CG2O_PARDISO_LIBRARY
  NAMES
    solver_eigen_pardiso
    cg2o_solver_eigen_pardiso
    cg2o_solver_pardiso
  HINTS ${_CG2O_HINTS}
  PATH_SUFFIXES lib lib64
)

# ------------------------------------------------------------
# Optional version extraction
# ------------------------------------------------------------

set(CG2O_VERSION "")
if(CG2O_INCLUDE_DIR AND EXISTS "${CG2O_INCLUDE_DIR}/cg2o/version.h")
  file(STRINGS "${CG2O_INCLUDE_DIR}/cg2o/version.h" _cg2o_version_major
       REGEX "^#define[ \t]+CG2O_VERSION_MAJOR[ \t]+[0-9]+")
  file(STRINGS "${CG2O_INCLUDE_DIR}/cg2o/version.h" _cg2o_version_minor
       REGEX "^#define[ \t]+CG2O_VERSION_MINOR[ \t]+[0-9]+")
  file(STRINGS "${CG2O_INCLUDE_DIR}/cg2o/version.h" _cg2o_version_patch
       REGEX "^#define[ \t]+CG2O_VERSION_PATCH[ \t]+[0-9]+")

  string(REGEX REPLACE ".*CG2O_VERSION_MAJOR[ \t]+([0-9]+).*" "\\1"
         _cg2o_major "${_cg2o_version_major}")
  string(REGEX REPLACE ".*CG2O_VERSION_MINOR[ \t]+([0-9]+).*" "\\1"
         _cg2o_minor "${_cg2o_version_minor}")
  string(REGEX REPLACE ".*CG2O_VERSION_PATCH[ \t]+([0-9]+).*" "\\1"
         _cg2o_patch "${_cg2o_version_patch}")

  if(_cg2o_major MATCHES "^[0-9]+$")
    set(CG2O_VERSION "${_cg2o_major}.${_cg2o_minor}.${_cg2o_patch}")
  endif()
endif()

# ------------------------------------------------------------
# Component handling
# ------------------------------------------------------------

if(NOT CG2O_FIND_COMPONENTS)
  set(CG2O_FIND_COMPONENTS core)
endif()

set(_CG2O_REQUIRED_VARS CG2O_INCLUDE_DIR)
set(_CG2O_ALL_LIBS)

foreach(_comp IN LISTS CG2O_FIND_COMPONENTS)
  if(_comp STREQUAL "core")
    list(APPEND _CG2O_REQUIRED_VARS CG2O_CORE_LIBRARY)
    if(CG2O_CORE_LIBRARY)
      list(APPEND _CG2O_ALL_LIBS "${CG2O_CORE_LIBRARY}")
    endif()

  elseif(_comp STREQUAL "cg2o")
    list(APPEND _CG2O_REQUIRED_VARS CG2O_CG2O_LIBRARY)
    if(CG2O_CG2O_LIBRARY)
      list(APPEND _CG2O_ALL_LIBS "${CG2O_CG2O_LIBRARY}")
    endif()

  elseif(_comp STREQUAL "umfpack")
    list(APPEND _CG2O_REQUIRED_VARS CG2O_UMFPACK_LIBRARY)
    if(CG2O_UMFPACK_LIBRARY)
      list(APPEND _CG2O_ALL_LIBS "${CG2O_UMFPACK_LIBRARY}")
    endif()

  elseif(_comp STREQUAL "pardiso")
    list(APPEND _CG2O_REQUIRED_VARS CG2O_PARDISO_LIBRARY)
    if(CG2O_PARDISO_LIBRARY)
      list(APPEND _CG2O_ALL_LIBS "${CG2O_PARDISO_LIBRARY}")
    endif()

  else()
    message(FATAL_ERROR
      "FindCG2O.cmake: unsupported component '${_comp}'. "
      "Supported components are: core, cg2o, umfpack, pardiso")
  endif()
endforeach()

list(REMOVE_DUPLICATES _CG2O_REQUIRED_VARS)
list(REMOVE_DUPLICATES _CG2O_ALL_LIBS)

find_package_handle_standard_args(CG2O
  REQUIRED_VARS ${_CG2O_REQUIRED_VARS}
  VERSION_VAR CG2O_VERSION
  HANDLE_COMPONENTS
)

# ------------------------------------------------------------
# Legacy-style variables
# ------------------------------------------------------------

if(CG2O_FOUND)
  set(CG2O_INCLUDE_DIRS "${CG2O_INCLUDE_DIR}")
  set(CG2O_LIBRARIES "${_CG2O_ALL_LIBS}")
endif()

# ------------------------------------------------------------
# Imported targets
# ------------------------------------------------------------

if(CG2O_FOUND AND NOT TARGET CG2O::core AND CG2O_CORE_LIBRARY)
  add_library(CG2O::core UNKNOWN IMPORTED)
  set_target_properties(CG2O::core PROPERTIES
    IMPORTED_LOCATION "${CG2O_CORE_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${CG2O_INCLUDE_DIR}"
  )
endif()

if(CG2O_FOUND AND NOT TARGET CG2O::umfpack AND CG2O_UMFPACK_LIBRARY)
  add_library(CG2O::umfpack UNKNOWN IMPORTED)
  set_target_properties(CG2O::umfpack PROPERTIES
    IMPORTED_LOCATION "${CG2O_UMFPACK_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${CG2O_INCLUDE_DIR}"
  )
  if(TARGET CG2O::core)
    set_property(TARGET CG2O::umfpack APPEND PROPERTY
      INTERFACE_LINK_LIBRARIES CG2O::core)
  endif()
endif()

if(CG2O_FOUND AND NOT TARGET CG2O::pardiso AND CG2O_PARDISO_LIBRARY)
  add_library(CG2O::pardiso UNKNOWN IMPORTED)
  set_target_properties(CG2O::pardiso PROPERTIES
    IMPORTED_LOCATION "${CG2O_PARDISO_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${CG2O_INCLUDE_DIR}"
  )
  if(TARGET CG2O::core)
    set_property(TARGET CG2O::pardiso APPEND PROPERTY
      INTERFACE_LINK_LIBRARIES CG2O::core)
  endif()
endif()

if(CG2O_FOUND AND NOT TARGET CG2O::cg2o AND CG2O_CG2O_LIBRARY)
  add_library(CG2O::cg2o UNKNOWN IMPORTED)
  set_target_properties(CG2O::cg2o PROPERTIES
    IMPORTED_LOCATION "${CG2O_CG2O_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${CG2O_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(
  CG2O_INCLUDE_DIR
  CG2O_CORE_LIBRARY
  CG2O_CG2O_LIBRARY
  CG2O_UMFPACK_LIBRARY
  CG2O_PARDISO_LIBRARY
)
