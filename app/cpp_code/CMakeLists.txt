cmake_minimum_required(VERSION 3.22)

project(acc_control_mpc_cg2o)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(cmake/CompilerSetting.cmake)

find_package(fmt REQUIRED)
find_package(g2o REQUIRED)
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

find_package(cg2o REQUIRED)

if(TARGET cg2o::cg2o)
    message(STATUS "Found cg2o main library")
else()
    message(FATAL_ERROR "Could not find cg2o main library")
endif()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../lib
                 ${CMAKE_CURRENT_BINARY_DIR}/lib_build)

set(EXC "mpc_g2o")
set(EXEC_NAME "${EXC}_exc")

# Default optimizer, overridable from outside
# using a cmake argument -DMPC_CG2O_OPTIMIZER=AL|ISPD|BIPM:
 set(MPC_CG2O_OPTIMIZER "BIPM" CACHE STRING "Optimizer backend for mpc executable")
set_property(CACHE MPC_CG2O_OPTIMIZER PROPERTY STRINGS AL ISPD BIPM)

function(set_mpc_optimizer target optimizer_var)
    if(${optimizer_var} STREQUAL "AL")
        target_compile_definitions(${target} PRIVATE USE_AL)
    elseif(${optimizer_var} STREQUAL "ISPD")
        target_compile_definitions(${target} PRIVATE USE_ISPD)
    elseif(${optimizer_var} STREQUAL "BIPM")
        target_compile_definitions(${target} PRIVATE USE_BIPM)
    else()
        message(FATAL_ERROR
            "Unknown optimizer '${${optimizer_var}}' for target ${target}. "
            "Use AL, ISPD, or BIPM")
    endif()
endfunction()

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
    ${Python3_INCLUDE_DIRS}
)

# Local optimizer selection for this executable only
set_mpc_optimizer(${EXEC_NAME} MPC_CG2O_OPTIMIZER)

if(MPC_USE_EXACT_JACOBIANS)
    target_compile_definitions(${EXEC_NAME} PRIVATE MPC_USE_EXACT_JACOBIANS)
    message(STATUS "Using Exact Jacobian")
endif()

if(TARGET cg2o::solver_eigen_umfpack)
    target_compile_definitions(${EXEC_NAME} PRIVATE MPC_USE_SOLVERS_UMFPACK)
    message(STATUS "Found cg2o UMFPACK solver")
else()
    message(STATUS "UMFPACK solver not found, will not be used")
endif()

if(TARGET cg2o::solver_eigen_pardiso)
    target_compile_definitions(${EXEC_NAME} PRIVATE MPC_USE_SOLVERS_PARDISO)
    target_link_libraries(${EXEC_NAME} PRIVATE cg2o::solver_eigen_pardiso)
    message(STATUS "Found cg2o PARDISO solver")
endif()

 

if(TARGET cg2o::solvers_eigen)
    target_compile_definitions(${EXEC_NAME} PRIVATE MPC_USE_SOLVERS_EIGEN)
    target_link_libraries(${EXEC_NAME} PRIVATE cg2o::solvers_eigen)
    message(STATUS "Found cg2o Eigen solvers")
endif()


if(TARGET g2o::solver_cholmod AND TARGET g2o::solver_csparse)
    target_compile_definitions(${EXEC_NAME} PRIVATE MPC_USE_G2O_SOLVERS)
    target_link_libraries(${EXEC_NAME} PRIVATE g2o::solver_cholmod  g2o::solver_csparse)
    message(STATUS "Found cg2o Eigen solvers")
endif()

 

target_compile_definitions(${EXEC_NAME} PRIVATE MPC_FEASIBLE_INITIALIZATION)

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
