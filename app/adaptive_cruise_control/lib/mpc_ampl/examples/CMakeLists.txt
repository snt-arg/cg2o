# as independent project
cmake_minimum_required(VERSION 3.13)
project(AmplCppExample)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Path to AMPL API root
set(AMPLAPI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../../../../../../dep/ampl.linux)
 
# Include paths
include_directories(${AMPLAPI_DIR}/amplapi/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/mpc_problem)

# Link directories
link_directories(${AMPLAPI_DIR}/amplapi/lib)

# Add executable
add_executable(ampl_example1 example1.cpp)

# Define compile-time constants for this target
target_compile_definitions(ampl_example1 PUBLIC
    AMPLAPI_DIR="${AMPLAPI_DIR}"
    MPC_AMPL_PATH="${CMAKE_CURRENT_SOURCE_DIR}"
)

# Link AMPL library
target_link_libraries(ampl_example1 ampl)

