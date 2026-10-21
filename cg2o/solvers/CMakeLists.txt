add_subdirectory(umfpack)

# Quasi-definite LDL^T solver (Eigen SimplicialLDLT). Always built,
# independent of CG2O_BUILD_EIGEN_SOLVERS: it only needs Eigen + g2o. The
# header can be used on its own; the library adds the g2o factory
# registration (gn_var_qdldlt, ...).
add_library(solver_qdldlt ${G2O_LIB_TYPE}
  eigen/solver_eigen_qdldlt.cpp
  eigen/linear_solver_eigen_qdldlt.h
)
target_include_directories(solver_qdldlt PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(solver_qdldlt PUBLIC
  g2o::core
  g2o::stuff
  ${CG2O_EIGEN3_EIGEN_TARGET}
)
set_target_properties(solver_qdldlt PROPERTIES
  OUTPUT_NAME ${CG2O_LIB_PREFIX}solver_qdldlt
  VERSION ${CG2O_LIB_VERSION}
  SOVERSION ${CG2O_LIB_SOVERSION}
)
target_compile_features(solver_qdldlt PUBLIC cxx_std_20)
install(TARGETS solver_qdldlt
  EXPORT ${CG2O_TARGETS_EXPORT_NAME}
  RUNTIME DESTINATION ${RUNTIME_DESTINATION}
  LIBRARY DESTINATION ${LIBRARY_DESTINATION}
  ARCHIVE DESTINATION ${ARCHIVE_DESTINATION}
  INCLUDES DESTINATION ${INCLUDES_DESTINATION}
)
if(NOT CG2O_BUILD_EIGEN_SOLVERS)  # otherwise eigen/ installs all its headers
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/eigen/linear_solver_eigen_qdldlt.h
          DESTINATION ${INCLUDES_INSTALL_DIR}/solvers/eigen)
endif()

# CHOLMOD with cached, pre-permuted analysis (built if CHOLMOD is found)
add_subdirectory(cholmod)

if(CG2O_BUILD_EIGEN_SOLVERS)
   add_subdirectory(eigen)
 endif()
 
if(CG2O_BUILD_PARDISO)
   add_subdirectory(pardiso)
endif()

