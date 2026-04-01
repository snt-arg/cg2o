set(SOURCE_FILES 
  solver_eigen_dense_full_lu.cpp 
  solver_eigen_dense_partial_lu.cpp
  solver_eigen_dense_qr.cpp
  solver_eigen_svd.cpp 
  solver_eigen_lu.cpp 
  solver_eigen_dense_qr.cpp 
  solver_eigen_lu.cpp
  solver_eigen_bicgstab.cpp
  solver_eigen_lscg.cpp
  solver_eigen_ldlt.cpp  
)

set(HEADER_FILES
  linear_solver_eigen_lu.h 
  linear_solver_eigen_qr.h 
  linear_solver_eigen_ldlt.h   
  linear_solver_eigen_dense_full_lu.h 
  linear_solver_eigen_dense_partial_lu.h
  linear_solver_eigen_dense_qr.h
  linear_solver_eigen_svd.h 
  linear_solver_eigen_bicgstab.h
  linear_solver_eigen_lscg.h
)

add_library(solvers_eigen ${G2O_LIB_TYPE}
 ${SOURCE_FILES}
 ${HEADER_FILES}
)

set_target_properties(solvers_eigen PROPERTIES OUTPUT_NAME ${LIB_PREFIX}solvers_eigen)
set_target_properties(solvers_eigen PROPERTIES
  VERSION ${CG2O_LIB_VERSION}
  SOVERSION ${CG2O_LIB_SOVERSION})
if (APPLE)
  set_target_properties(solvers_eigen PROPERTIES INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif()
target_compile_features(solvers_eigen PUBLIC cxx_std_17)
target_link_libraries(solvers_eigen PUBLIC ${G2O_LIBRARIES} ${CG2O_EIGEN3_EIGEN_TARGET})

target_include_directories(solvers_eigen PUBLIC
  "$<BUILD_INTERFACE:${cg2o_SOURCE_DIR};${PROJECT_BINARY_DIR}>"
  $<INSTALL_INTERFACE:include/cg2o/solvers/eigen>
)

install(TARGETS solvers_eigen
  EXPORT ${CG2O_TARGETS_EXPORT_NAME}
  RUNTIME DESTINATION ${RUNTIME_DESTINATION}
  LIBRARY DESTINATION ${LIBRARY_DESTINATION}
  ARCHIVE DESTINATION ${ARCHIVE_DESTINATION}
  INCLUDES DESTINATION ${INCLUDES_DESTINATION}
)

file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
install(FILES ${headers} DESTINATION ${INCLUDES_INSTALL_DIR}/solvers/eigen)
