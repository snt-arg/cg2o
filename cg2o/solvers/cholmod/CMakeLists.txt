# LinearSolverCholmodFast: CHOLMOD (SuiteSparse) with the ordering and
# symbolic analysis kept across optimize() calls and a pre-permuted matrix.
find_path(CHOLMOD_INCLUDE_DIR
  NAMES cholmod.h
  PATHS /usr/include/suitesparse /usr/include
)
find_library(CHOLMOD_LIB
  NAMES cholmod
)
if(NOT CHOLMOD_INCLUDE_DIR OR NOT CHOLMOD_LIB)
  message(STATUS "CHOLMOD not found (cholmod.h / libcholmod): LinearSolverCholmodFast is not built")
  return()
endif()

set(SOURCE_FILES
  solver_cholmod_fast.cpp
)

set(HEADER_FILES
  linear_solver_cholmod_fast.h
)

add_library(solver_cholmod_fast ${G2O_LIB_TYPE}
  ${SOURCE_FILES}
  ${HEADER_FILES}
)

target_include_directories(solver_cholmod_fast
  PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    ${CHOLMOD_INCLUDE_DIR}
)

target_link_libraries(solver_cholmod_fast
  PUBLIC
    ${CG2O_EIGEN3_EIGEN_TARGET}
    g2o::core
    g2o::stuff
    ${CHOLMOD_LIB}
)

set_target_properties(solver_cholmod_fast PROPERTIES
  OUTPUT_NAME ${CG2O_LIB_PREFIX}solver_cholmod_fast
  VERSION ${CG2O_LIB_VERSION}
  SOVERSION ${CG2O_LIB_SOVERSION}
)

if(APPLE)
  set_target_properties(solver_cholmod_fast PROPERTIES INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif()

target_compile_features(solver_cholmod_fast PUBLIC cxx_std_20)

install(TARGETS solver_cholmod_fast
  EXPORT ${CG2O_TARGETS_EXPORT_NAME}
  RUNTIME DESTINATION ${RUNTIME_DESTINATION}
  LIBRARY DESTINATION ${LIBRARY_DESTINATION}
  ARCHIVE DESTINATION ${ARCHIVE_DESTINATION}
  INCLUDES DESTINATION ${INCLUDES_DESTINATION}
)

file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
install(FILES ${headers} DESTINATION ${INCLUDES_INSTALL_DIR}/solvers/cholmod)
