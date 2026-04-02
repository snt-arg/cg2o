set(source
  sparse_optimizer_bipm.cpp
  sparse_optimizer_ispd.cpp
  sparse_optimizer_al.cpp
)

set(header
  base_fixed_sized_edge_eq.h
  base_fixed_sized_edge_ineq.h
  sparse_optimizer_al.h
  sparse_optimizer_bipm.h 
  sparse_optimizer_ispd.h
  sparse_optimizer_for_constraints.h
  sparse_optimizer_for_constraints.hpp
  termination_criteria.h
  vertex_lagrange_multiplier.h
)

# Define the library
add_library(core ${G2O_LIB_TYPE} ${source} ${header})

# Include directories
target_include_directories(core PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_link_libraries(core PUBLIC
   g2o::core ${CG2O_EIGEN3_EIGEN_TARGET}
   g2o::stuff
)

set_target_properties(core PROPERTIES 
  OUTPUT_NAME ${CG2O_LIB_PREFIX}core 
  VERSION ${CG2O_LIB_VERSION}
  SOVERSION ${CG2O_LIB_SOVERSION}
)

target_compile_features(core PUBLIC cxx_std_20)

if (APPLE)
  set_target_properties(core PROPERTIES INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif()

install(TARGETS core
  EXPORT ${CG2O_TARGETS_EXPORT_NAME}
  RUNTIME DESTINATION ${RUNTIME_DESTINATION}
  LIBRARY DESTINATION ${LIBRARY_DESTINATION}
  ARCHIVE DESTINATION ${ARCHIVE_DESTINATION}
  INCLUDES DESTINATION ${INCLUDES_DESTINATION}
)

file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

install(FILES ${headers} DESTINATION ${INCLUDES_INSTALL_DIR}/core)
