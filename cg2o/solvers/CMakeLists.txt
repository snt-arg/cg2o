add_subdirectory(umfpack)


if(CG2O_BUILD_EIGEN_SOLVERS)
   add_subdirectory(eigen)
 endif()
 
if(CG2O_BUILD_PARDISO)
   add_subdirectory(pardiso)
endif()

