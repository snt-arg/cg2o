cmake_minimum_required(VERSION 3.22)

project(acc_control_mpc_cg2o)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(cmake/CompilerSetting.cmake) # Compiler settings

#list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(fmt REQUIRED)
find_package(g2o  REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    add_definitions(-DG2O_USE_OPENMP)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")


find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
include_directories(${Python3_INCLUDE_DIRS})
 

#set(CG2O_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../build")
find_package(cg2o REQUIRED)

if (TARGET cg2o::cg2o)
    message(STATUS "Found cg2o main library")
else()
    message(FATAL_ERROR "Could not find cg2o  main library")
endif()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../lib
                 ${CMAKE_CURRENT_BINARY_DIR}/lib_build) # here the controller is built as a library and linked to the executable

set(EXC "mpc_space")
set(EXEC_NAME "${EXC}_exc")  # mpc_space_exc  
 
add_executable(${EXEC_NAME} src/main.cpp)

target_link_libraries(${EXEC_NAME} PRIVATE
  g2o::core
  g2o::stuff
  Eigen3::Eigen
  fmt::fmt
  cg2o::cg2o
  mpc_cg2o_lib
  OpenMP::OpenMP_CXX
  ${Python3_LIBRARIES}
  yaml-cpp
)

target_include_directories(${EXEC_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../lib/data_saver
    ${CMAKE_CURRENT_SOURCE_DIR}/../lib
)

target_compile_definitions(${EXEC_NAME} PRIVATE USE_ISPD)




 
if(MPC_USE_EXACT_JACOBIANS)
target_compile_definitions(${EXEC_NAME} PRIVATE MPC_USE_EXACT_JACOBIANS)
message(STATUS  "Using Exact Jacobian")
endif()
 
 



if(TARGET cg2o::solver_eigen_umfpack)
     add_compile_definitions(USE_SOLVERS_UMFPACK) 
    # no need to link the UMFPACK solver directly, as it's part of the cg2o target
    message(STATUS "Found cg2o UMFPACK solver")
    else()  
    message(STATUS "UMFPACK solver not found, will not be used")
endif()

if(TARGET cg2o::solver_eigen_pardiso)
    add_compile_definitions(USE_SOLVERS_PARDISO)
    target_link_libraries(${EXEC_NAME} PRIVATE cg2o::solver_eigen_pardiso)
    message(STATUS "Found cg2o PARDISO solver")
endif()

if(TARGET cg2o::solvers_eigen)
    add_compile_definitions(USE_SOLVERS_EIGEN)
    target_link_libraries(${EXEC_NAME} PRIVATE cg2o::solver_eigen_pardiso)
    message(STATUS "Found cg2o Eingen solvers")
endif()
   






 
 

# save the executable to the bin 
set_target_properties(${EXEC_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/results/data"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_CURRENT_SOURCE_DIR}/results/data"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_CURRENT_SOURCE_DIR}/results/data"
    RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL "${CMAKE_CURRENT_SOURCE_DIR}/results/data"
    RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO "${CMAKE_CURRENT_SOURCE_DIR}/results/data"
)

 
set(DEFAULT_CONFIG_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/config/config.yaml"
)

target_compile_definitions(${EXEC_NAME} PRIVATE
  DEFAULT_CONFIG_PATH="${DEFAULT_CONFIG_PATH}"
)

target_compile_definitions(${EXEC_NAME} PRIVATE
  DEFAULT_CONFIG_PATH="${DEFAULT_CONFIG_PATH}"
)


