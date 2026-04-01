

#  min ||xy(1)-a||^2 +  ||xy(2)-b||^2 + ||z-c||^2 
#  s.t.    xy(2) + z == 3, xy(1)+ xy(2) == 2.5      xy <= {1, 4}   z <= 5
# add executable test
add_executable(example_1 src/example_1.cpp)
target_link_libraries(example_1 ${G2O_LIBRARIES} cg2o Eigen3::Eigen fmt::fmt)
target_include_directories(example_1 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if(CG2O_BUILD_EIGEN_SOLVERS)
 target_link_libraries(example_1 solvers_eigen_extension)
 endif()
 
if(CG2O_BUILD_PARDISO)
 target_link_libraries(example_1 solver_eigen_pardiso)
endif()





# add executable test
# minimize     0.5 * ||x1-1||^2 + 0.5 * ||x2-1||^2
#    subject to    x1^2 + x2^2 - 1 = 0
#                  4*x1^2 + 0.25*x2^2 - 1 <= 0
add_executable(example_2 src/example_2.cpp)
target_link_libraries(example_2 ${G2O_LIBRARIES} cg2o Eigen3::Eigen fmt::fmt)
target_include_directories(example_2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)


#     minimize     0.5*||x1 + e^(-x2)||^2 + 0.5*||x1^2 + 2*x2 + 1||^2
#    subject to    x1 + x1^3 + x2 + x2^2  = 0                
# add executable test
add_executable(example_3 src/example_3.cpp)
target_link_libraries(example_3 ${G2O_LIBRARIES} cg2o Eigen3::Eigen fmt::fmt)
target_include_directories(example_3 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
