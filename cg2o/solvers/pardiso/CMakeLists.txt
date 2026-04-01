set(SOURCE_FILES 
  solver_eigen_pardiso_lu.cpp
)

set(HEADER_FILES
  linear_solver_eigen_pardiso_lu.h
)

add_library(solver_eigen_pardiso ${G2O_LIB_TYPE}
  ${SOURCE_FILES}
  ${HEADER_FILES}
)

target_include_directories(solver_eigen_pardiso PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
  
  
  target_include_directories(solver_eigen_pardiso PRIVATE  
  ${MKL_INCLUDE_DIR} 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  g2o::stuff
  g2o::core 
  )
  
  
target_link_libraries(solver_eigen_pardiso PRIVATE 
  g2o::stuff
  g2o::core   
  ${CG2O_EIGEN3_EIGEN_TARGET}
)

target_link_libraries(solver_eigen_pardiso PUBLIC 
MKL 
)

# Define preprocessor macros
target_compile_definitions(solver_eigen_pardiso PRIVATE
EIGEN_USE_MKL_ALL
)


set_target_properties(solver_eigen_pardiso PROPERTIES
BUILD_RPATH "${MKL_LIBRARY_DIR}"
)


set_target_properties(solver_eigen_pardiso PROPERTIES OUTPUT_NAME ${CG2O_LIB_PREFIX}solver_eigen_pardiso)
set_target_properties(solver_eigen_pardiso PROPERTIES
VERSION ${CG2O_LIB_VERSION}
SOVERSION ${CG2O_LIB_SOVERSION})
if (APPLE)
  set_target_properties(solver_eigen_pardiso PROPERTIES INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif()
target_compile_features(solver_eigen_pardiso PUBLIC cxx_std_17)



install(TARGETS solver_eigen_pardiso
  EXPORT ${CG2O_TARGETS_EXPORT_NAME}
  RUNTIME DESTINATION ${RUNTIME_DESTINATION}
  LIBRARY DESTINATION ${LIBRARY_DESTINATION}
  ARCHIVE DESTINATION ${ARCHIVE_DESTINATION}
  INCLUDES DESTINATION ${INCLUDES_DESTINATION}
)

file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
install(FILES ${headers} DESTINATION ${INCLUDES_INSTALL_DIR}/solvers/pardiso)
 
 