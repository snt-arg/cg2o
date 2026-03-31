add_subdirectory(solvers)
add_subdirectory(core)

add_library(cg2o INTERFACE)

target_link_libraries(cg2o
  INTERFACE
  core
  solver_eigen_umfpack
)
 
 
 
if(CG2O_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
