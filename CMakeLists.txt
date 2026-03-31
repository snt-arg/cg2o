cmake_minimum_required(VERSION 3.22)
project(cg2o LANGUAGES CXX)

set(CMAKE_LEGACY_CYGWIN_WIN32 0)
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15" CACHE STRING "Minimum OS X deployment target")

# The library prefix
set(LIB_PREFIX cg2o_)

set(g2o_C_FLAGS)
set(g2o_CXX_FLAGS)

set(CG2O_LIB_VERSION   "0.1.0" CACHE STRING "cg2o library version")
set(CG2O_LIB_SOVERSION "0.1"   CACHE STRING "cg2o library soversion")
set(CG2O_VERSION 1.0.0)

#----------------------------------------------------------------------------------------
##options 
# Configuration option: Default to ISPD optimizer, but allow switching to AL or IPM
option(CG2O_BUILD_EXAMPLES "Build the examples" ON)
option(CG2O_BUILD_MPC "Build the MPC controller" OFF)
option(USE_G2O_SOLVERS "Enable g2o-based solvers (links g2o + csparse and cholmod)" OFF)
option(CG2O_BUILD_PARDISO "allow using the Paradiso solver" ON)
option(CG2O_BUILD_EIGEN_SOLVERS "allow building the Eigen solver" OFF)

set(CG2O_OPTIMIZER "BIPM" CACHE STRING "Optimizer type: ISPD, AL,  or BIPM")
set_property(CACHE CG2O_OPTIMIZER PROPERTY STRINGS ISPD AL BIPM)




#----------------------------------------------------------------------------------------
# Find required packages
find_package(Eigen3 3.3 REQUIRED)
if (TARGET Eigen3::Eigen)
  set(CG2O_EIGEN3_EIGEN_TARGET Eigen3::Eigen)
else()
  include_directories(${EIGEN3_INCLUDE_DIR})
endif ()
find_package(fmt REQUIRED)
find_package(g2o 1.0.0 REQUIRED)




# Specify the path to g2o include directory and libraries
set(G2O_INCLUDE_DIRS "/usr/include/include" CACHE PATH "g2o include directory")
set(G2O_LIBRARY_DIRS "/usr/include/lib" CACHE PATH "g2o library directory")
set(G2O_LIBRARIES
    g2o::core
    g2o::stuff
)

if(USE_G2O_SOLVERS)
    add_compile_definitions(USE_G2O_SOLVERS)

    list(APPEND G2O_LIBRARIES
        g2o::solver_csparse
        g2o::csparse_extension
        g2o::solver_cholmod
    )
endif()


# Manually specify g2o include directories and libraries
include_directories(${G2O_INCLUDE_DIRS})
link_directories(${G2O_LIBRARY_DIRS})


# OpenMP support
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    add_compile_definitions(G2O_USE_OPENMP)
endif()


if(CG2O_BUILD_PARDISO)
    add_compile_definitions(CG2O_BUILD_PARDISO)
    # Set MKL paths
    set(MKL_ROOT "/opt/intel/oneapi/mkl/latest")
    set(MKL_INCLUDE_DIR "${MKL_ROOT}/include")
    set(MKL_LIBRARY_DIR "${MKL_ROOT}/lib/intel64")

    # Option to choose threading model (TBB or OpenMP)
    option(MKL_USE_TBB "Use MKL with TBB threading" OFF)
    if(MKL_USE_TBB)
        find_package(TBB REQUIRED)
        set(MKL_THREAD_LIB "-lmkl_tbb_thread")
        set(THREAD_LIB "TBB::tbb")
        message(STATUS "Using MKL with TBB threading for Pardiso solver")
    else()
        set(MKL_THREAD_LIB "-lmkl_gnu_thread")
        set(THREAD_LIB "-lgomp")
        message(STATUS "Using MKL with OpenMP threading for Pardiso solver")
    endif()

    # Define MKL library
    add_library(MKL INTERFACE)
    target_include_directories(MKL INTERFACE ${MKL_INCLUDE_DIR})
    target_link_libraries(MKL INTERFACE
        -L${MKL_LIBRARY_DIR}
        -lmkl_intel_lp64
        ${MKL_THREAD_LIB}
        -lmkl_core
        -lmkl_rt
        ${THREAD_LIB}
        -lpthread -lm -ldl
    )
    set(CMAKE_BUILD_RPATH "${MKL_LIBRARY_DIR}")
endif()



# Manually specify bipm_g2o include directories and libraries
include_directories(${G2O_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include)  
link_directories(${G2O_LIBRARY_DIRS})


if(CG2O_BUILD_EXAMPLES)
    message(STATUS "Building with examples enabled")
    add_compile_definitions(CG2O_BUILD_EXAMPLES)
endif()


if(CG2O_OPTIMIZER STREQUAL "AL")
  add_compile_definitions(USE_AL)
elseif(CG2O_OPTIMIZER STREQUAL "ISPD")
  add_compile_definitions(USE_ISPD)
elseif(CG2O_OPTIMIZER STREQUAL "BIPM")
  add_compile_definitions(USE_BIPM)
else()
  message(FATAL_ERROR "Unknown CG2O_OPTIMIZER=${CG2O_OPTIMIZER}. Use AL, ISPD, or BIPM")
endif()


# add the cg2o extenstion
add_subdirectory(cg2o)

if(CG2O_BUILD_MPC)

    option(MPC_USE_EXACT_JACOBIANS "Use analytical Jacobians instead of numerical ones" ON)
    option(MPC_AMPL_LANGUAGE "Enable AMPL-based MPC formulation instead of G2O"   OFF)
    option(MPC_BOTH_MODELS "Enable G2O"  OFF)

    
    if(MPC_USE_EXACT_JACOBIANS)
        add_compile_definitions(MPC_USE_EXACT_JACOBIANS)
        message(STATUS  "Using analytical Jacobian")
    endif()

    if(MPC_BOTH_MODELS)
        add_compile_definitions(MPC_BOTH_MODELS)
        message(STATUS  "Using both models")
    endif()
    
    # either MPC_AMPL_LANGUAGE or USE_CASADI_LANGUAGE or none of them
    if(MPC_AMPL_LANGUAGE AND MPC_CASADI_TOOL)
        message(FATAL_ERROR "Cannot enable both MPC_AMPL_LANGUAGE and MPC_CASADI_TOOL at the same time.")
    endif()


    add_subdirectory(mpc)  

endif()


