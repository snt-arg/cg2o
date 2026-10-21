add_subdirectory(umfpack)

# Header-only quasi-definite LDL^T solver (Eigen SimplicialLDLT). Always
# available, independent of CG2O_BUILD_EIGEN_SOLVERS: it only needs Eigen + g2o.
add_library(solver_qdldlt INTERFACE)
target_include_directories(solver_qdldlt INTERFACE
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(solver_qdldlt INTERFACE
  g2o::core
  g2o::stuff
  ${CG2O_EIGEN3_EIGEN_TARGET}
)
target_compile_features(solver_qdldlt INTERFACE cxx_std_20)
install(TARGETS solver_qdldlt EXPORT ${CG2O_TARGETS_EXPORT_NAME})
if(NOT CG2O_BUILD_EIGEN_SOLVERS)  # otherwise eigen/ installs all its headers
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/eigen/linear_solver_eigen_qdldlt.h
          DESTINATION ${INCLUDES_INSTALL_DIR}/solvers/eigen)
endif()


if(CG2O_BUILD_EIGEN_SOLVERS)
   add_subdirectory(eigen)
 endif()
 
if(CG2O_BUILD_PARDISO)
   add_subdirectory(pardiso)
endif()

