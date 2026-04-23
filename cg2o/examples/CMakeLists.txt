# Per-example optimizer selection
set(EXAMPLE_1_OPTIMIZER "ISPD"   CACHE STRING "Optimizer for example_1")
set(EXAMPLE_2_OPTIMIZER "ISPD"   CACHE STRING "Optimizer for example_2")
set(EXAMPLE_3_OPTIMIZER "ISPD"   CACHE STRING "Optimizer for example_3")

set_property(CACHE EXAMPLE_1_OPTIMIZER PROPERTY STRINGS AL ISPD BIPM)
set_property(CACHE EXAMPLE_2_OPTIMIZER PROPERTY STRINGS AL ISPD BIPM)
set_property(CACHE EXAMPLE_3_OPTIMIZER PROPERTY STRINGS AL ISPD BIPM)

function(set_example_optimizer target optimizer)
  if(${optimizer} STREQUAL "AL")
    target_compile_definitions(${target} PRIVATE USE_AL)
  elseif(${optimizer} STREQUAL "ISPD")
    target_compile_definitions(${target} PRIVATE USE_ISPD)
  elseif(${optimizer} STREQUAL "BIPM")
    target_compile_definitions(${target} PRIVATE USE_BIPM)
  else()
    message(FATAL_ERROR "Unknown optimizer '${${optimizer}}' for target ${target}. Use AL, ISPD, or BIPM")
  endif()
endfunction()


#  min ||xy(1)-a||^2 +  ||xy(2)-b||^2 + ||z-c||^2
#  s.t.    xy(2) + z == 3, xy(1)+ xy(2) == 2.5      xy <= {1, 4}   z <= 5
add_executable(example_1 src/example_1.cpp)
target_link_libraries(example_1 PRIVATE ${G2O_LIBRARIES} cg2o Eigen3::Eigen fmt::fmt)
target_include_directories(example_1 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_example_optimizer(example_1 EXAMPLE_1_OPTIMIZER)

if(CG2O_BUILD_EIGEN_SOLVERS)
  target_link_libraries(example_1 PRIVATE solvers_eigen)
endif()

if(CG2O_BUILD_PARDISO)
  target_link_libraries(example_1 PRIVATE solver_eigen_pardiso)
endif()


# minimize     0.5 * ||x1-1||^2 + 0.5 * ||x2-1||^2
# subject to   x1^2 + x2^2 - 1 = 0
#              4*x1^2 + 0.25*x2^2 - 1 <= 0
add_executable(example_2 src/example_2.cpp)
target_link_libraries(example_2 PRIVATE ${G2O_LIBRARIES} cg2o Eigen3::Eigen fmt::fmt)
target_include_directories(example_2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_example_optimizer(example_2 EXAMPLE_2_OPTIMIZER)


# minimize     0.5*||x1 + e^(-x2)||^2 + 0.5*||x1^2 + 2*x2 + 1||^2
# subject to   x1 + x1^3 + x2 + x2^2 = 0
add_executable(example_3 src/example_3.cpp)
target_link_libraries(example_3 PRIVATE ${G2O_LIBRARIES} cg2o Eigen3::Eigen fmt::fmt)
target_include_directories(example_3 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_example_optimizer(example_3 EXAMPLE_3_OPTIMIZER)
