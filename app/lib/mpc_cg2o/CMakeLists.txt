 

# Use it as library name
set(LIB_NAME "mpc_cg2o_lib")


# List of source files
file(GLOB mpc_edges_sources "${CMAKE_CURRENT_SOURCE_DIR}/mpc_edges/*.cpp")
file(GLOB mpc_edges_headers "${CMAKE_CURRENT_SOURCE_DIR}/mpc_edges/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/mpc_edges/*.hpp")

set(source 
)

# List of header files
set(header
  mpc_vertices.h
  mpc_vertices.hpp
  mpc_edges.h
  mpc_formulation.h
  mpc_parameters.h
 
)

# Define the library
add_library(${LIB_NAME} ${mpc_edges_sources}
  ${mpc_edges_headers}
  ${source}
  ${header}
)

 
# Link the library against its dependencies
target_include_directories(${LIB_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})  # Expose header files for dependent targets

target_compile_features(${LIB_NAME} PUBLIC cxx_std_17)


target_link_libraries(${LIB_NAME} PUBLIC ${G2O_LIBRARIES} Eigen3::Eigen cg2o::cg2o)


target_include_directories(${LIB_NAME}
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)



install(FILES
  ${header}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mpc_cg2o
)

install(DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/mpc_edges
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mpc_cg2o
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)

