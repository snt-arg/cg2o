cmake_minimum_required(VERSION 3.8)
project(acc_control_mpc_g2o)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()



# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(rclcpp_action REQUIRED)
find_package(acc_interfaces REQUIRED)
find_package(std_msgs REQUIRED)
find_package(yaml-cpp REQUIRED)



# Find required packages
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(fmt REQUIRED)
find_package(g2o 1.0.0 REQUIRED)
find_package(spdlog REQUIRED)
find_package(cg2o REQUIRED)


add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../../lib
                 ${CMAKE_CURRENT_BINARY_DIR}/lib_build) # here the controller is built as a library and linked to the executable
 


# OpenMP support
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    add_definitions(-DG2O_USE_OPENMP)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")


# Add the executable
add_executable(acc_control_mpc_g2o_ros src/acc_control_mpc_g2o_node.cpp)

  
# Add the include directory
target_include_directories(acc_control_mpc_g2o_ros PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)


# Link libraries
target_link_libraries(acc_control_mpc_g2o_ros
  g2o::core
  g2o::stuff
  Eigen3::Eigen
  fmt::fmt
  cg2o::cg2o
  mpc_cg2o_lib
  OpenMP::OpenMP_CXX
  yaml-cpp
 )

   

target_compile_definitions(acc_control_mpc_g2o_ros PRIVATE MPC_USE_EXACT_JACOBIANS)




# Add dependencies
ament_target_dependencies(acc_control_mpc_g2o_ros
  rclcpp
  rclcpp_action
  acc_interfaces
  yaml-cpp
)

# Install targets
install(TARGETS acc_control_mpc_g2o_ros
  DESTINATION lib/${PROJECT_NAME}
)

# Install include files
install(DIRECTORY include/
  DESTINATION include/
)

 

ament_package()
