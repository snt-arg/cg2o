add_subdirectory(core)
add_subdirectory(solvers)

add_library(cg2o INTERFACE)
target_link_libraries(cg2o INTERFACE core)

if(TARGET solver_eigen_umfpack)
  target_link_libraries(cg2o INTERFACE solver_eigen_umfpack)
endif()

add_library(cg2o::core ALIAS core)
add_library(cg2o::cg2o ALIAS cg2o)
if(TARGET solver_eigen_umfpack)
  add_library(cg2o::solver_eigen_umfpack ALIAS solver_eigen_umfpack)
endif()
if(TARGET solver_eigen_pardiso)
  add_library(cg2o::solver_eigen_pardiso ALIAS solver_eigen_pardiso)
endif()
if(TARGET solvers_eigen)
  add_library(cg2o::solvers_eigen ALIAS solvers_eigen)
endif()

install(TARGETS cg2o
  EXPORT ${CG2O_TARGETS_EXPORT_NAME}
)

if(CG2O_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
