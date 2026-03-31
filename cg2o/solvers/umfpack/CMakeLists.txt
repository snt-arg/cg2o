set(SOURCE_FILES
  solver_eigen_umfpack_lu.cpp
)

set(HEADER_FILES
  linear_solver_eigen_umfpack_lu.h
)

find_path(UMFPACK_INCLUDE_DIR
  NAMES umfpack.h
  PATHS /usr/include/suitesparse /usr/include
)

if(NOT UMFPACK_INCLUDE_DIR)
  message(FATAL_ERROR "UMFPACK headers not found (umfpack.h). Install libsuitesparse-dev.")
endif()

find_library(UMFPACK_LIB
  NAMES umfpack
)

if(NOT UMFPACK_LIB)
  message(FATAL_ERROR "UMFPACK library not found (libumfpack). Install libsuitesparse-dev.")
endif()

add_library(solver_eigen_umfpack ${G2O_LIB_TYPE}
  ${SOURCE_FILES}
  ${HEADER_FILES}
)

target_include_directories(solver_eigen_umfpack
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}   
    ${UMFPACK_INCLUDE_DIR}        
  PRIVATE
    ${G2O_INCLUDE_DIRS}           
)

 target_link_libraries(solver_eigen_umfpack PRIVATE
  ${G2O_LIBRARIES}
  ${CG2O_EIGEN3_EIGEN_TARGET}
  ${UMFPACK_LIB}
)


#----------------


set_target_properties(solver_eigen_umfpack PROPERTIES OUTPUT_NAME ${LIB_PREFIX}solver_eigen_umfpack)
set_target_properties(solver_eigen_umfpack PROPERTIES
  VERSION ${CG2O_LIB_VERSION}
  SOVERSION ${CG2O_LIB_SOVERSION})
if (APPLE)
  set_target_properties(solver_eigen_umfpack PROPERTIES INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif()
target_compile_features(solver_eigen_umfpack PUBLIC cxx_std_17)

target_include_directories(solver_eigen_umfpack PUBLIC
  "$<BUILD_INTERFACE:${cg2o_SOURCE_DIR};${PROJECT_BINARY_DIR}>"
  $<INSTALL_INTERFACE:include/cg2o/solvers/umfpack>
)

install(TARGETS solver_eigen_umfpack
  EXPORT ${CG2O_TARGETS_EXPORT_NAME}
  RUNTIME DESTINATION ${RUNTIME_DESTINATION}
  LIBRARY DESTINATION ${LIBRARY_DESTINATION}
  ARCHIVE DESTINATION ${ARCHIVE_DESTINATION}
  INCLUDES DESTINATION ${INCLUDES_DESTINATION}
)

file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
install(FILES ${headers} DESTINATION ${INCLUDES_INSTALL_DIR}/solvers/umfpack)

 



