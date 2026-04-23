#!/bin/bash


cat <<'EOF'

============================================================
CG2O Development Container
============================================================




## To build the CG2O examples, run:
====================================

  mkdir -p /ws/build && cd /ws/build && \
  rm -rf CMakeFiles  &&  rm -f CMakeCache.txt &&\
  cmake \
    -DCG2O_BUILD_EXAMPLES=ON \
    -DUSE_G2O_SOLVERS=ON \ 
    -DEXAMPLE_1_OPTIMIZER=ISPD \
    -DEXAMPLE_2_OPTIMIZER=AL \
    -DEXAMPLE_3_OPTIMIZER=BIPM \ 
    -DCG2O_BUILD_PARDISO=${ENABLE_PARDISO} \
    -DCG2O_BUILD_EIGEN_SOLVERS=ON \
    ..
  make -j$(nproc)

## Notes:
#  -DEXAMPLE_X_OPTIMIZER selects which optimizer is linked into the examples.
#  Options include: ISPD, BIPM, AL

## To run the examples:

  /ws/build/cg2o/examples/example_1 150 0 0 
  /ws/build/cg2o/examples/example_2
  /ws/build/cg2o/examples/example_3
  
  

## To build the MPC project:
====================================
  cd /ws/app/cpp_code  && \
  mkdir -p build && cd build && \
  rm -rf CMakeFiles  &&  rm -f CMakeCache.txt && \
  cd ../build 
  cmake .. \
    -DMPC_FEASIBLE_INITIALIZATION==OFF \
    -DMPC_CG2O_OPTIMIZER=ISPD\
    -DMPC_USE_NUMERICAL_JACOBIANS=OFF \
     && \
     make -j$(nproc)
  
   
 with   -DMPC_FEASIBLE_INITIALIZATION==ON \
     use       -DMPC_CG2O_OPTIMIZER=BIPM\


## To change the Controller and the Solver paramter modify the yaml on the path:
  /ws/app/cpp_code/config/config.yaml
## To run the MPC example:
  /ws/app/cpp_code/results/data/mpc_g2o_exc
 



## To build the ROS package:
===================================================
  source /opt/ros/jazzy/setup.bash  && \
  cd /ws/app/ros_ws  && \
  rm -rf build install log  && \
  colcon build  
  
## To run the ROS example:
 source  /ws/app/ros_ws/install/setup.bash


	ros2 run acc_control_mpc_g2o acc_control_mpc_g2o_ros --ros-args -p numberOfIterations:=40 -p horizon_length:=3 & \
	sleep 2; \
	ros2 action send_goal /acc_control_mpc_g2o_action acc_interfaces/action/AccControlMPC "{horizon_length: 3, v_p: 15.0, v_h: 12.0, a_p: 0.5, a_p_weighted: 0.5, d_h: 14.5, force_prev: 10.0}";\
	pkill -9 -f acc_control_mpc_g2o_ros

## Useful paths:
# Project root:   /ws
# C++ project:    /ws/app/cpp_code
# ROS workspace:  /ws/app/ros_ws
# ROS Eample:     /ws/app/ros_ws/sripts/run_ros_example.py

============================================================

EOF
    
