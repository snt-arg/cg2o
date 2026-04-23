# cg2o README.md

> **cg2o** — Constrained General Graph Optimization  
> A modern extension of **g2o** for **hard-constrained factor graph optimization**, introducing native support for **nonlinear inequality and equality constraints** using advanced optimization methods for robotics, control, and large-scale sparse systems.

---

## Overview

Factor graph optimization has become a standard tool in robotics for SLAM, state estimation, calibration, and trajectory optimization. However, classical frameworks such as g2o were primarily designed for unconstrained least-squares problems.

**cg2o** extends this paradigm to solve **general constrained nonlinear sparse optimization problems** while preserving the efficiency, modularity, and sparse structure of factor graphs.

The library introduces native constrained optimization solvers including:

- **ISPD-IPM** — Infeasible-Start Primal-Dual Interior Point Method (**main contribution of the paper**)
- **BIPM** — Barrier Interior Point Method
- **AL** — Augmented Lagrangian Method

This repository contains the implementation used in the paper:

**cg2o: A Constrained Factor Graph Solver via the Infeasible Start Interior Point Method for Nonlinear MPC**

---

## Why cg2o?

Most constrained optimization libraries operate on dense NLP formulations.

cg2o instead solves constrained problems directly on sparse factor graphs, enabling:

- Sparse scalable optimization
- Native graph abstraction
- Fast real-time MPC
- Equality + inequality constraints
- Reuse of existing g2o models
- Multiple constrained solvers in one framework

---

## Main Scientific Contribution

The key innovation is the implementation of the **ISPD-IPM constrained factor graph solver**.

### Traditional Barrier Methods Require

- Strictly feasible initialization
- Nested outer/inner loops
- Sensitive parameter tuning

### cg2o ISPD-IPM Provides

- Infeasible initialization
- Unified single optimization loop
- Dedicated inequality factors from KKT conditions
- Real-time performance

Reported in the paper: up to **58% fewer iterations** than a previous barrier-based factor graph solver.

---

## Available Solvers

```cpp
cg2o::SparseOptimizerISPD optimizer;   // Recommended
cg2o::SparseOptimizerBIPM optimizer;   // Barrier IPM
cg2o::SparseOptimizerAL optimizer;     // Augmented Lagrangian
```

---

## Supported Constraints

### Equality Constraints

```math
h(x)=0
```

Examples:

- system dynamics
- kinematic constraints
- calibration constraints

### Inequality Constraints

```math
g(x)\le 0
```

Examples:

- actuator limits
- collision avoidance
- speed bounds
- safety distances

---

## Applications

### Robotics

- SLAM with hard constraints
- Trajectory optimization
- Multi-robot systems
- Motion planning

### Control

- Nonlinear MPC
- Autonomous driving
- Drone control
- Energy management

---

## Build with Docker

### Build Image

```bash
docker build \
  --build-arg ENABLE_PARDISO=OFF \
  -t arg/cg2o \
  -f ../docker/Dockerfile ..
```

### Run Container

```bash
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

docker run -it --rm \
    --name cg2o \
    --user $(id -u):$(id -g) \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    -e DISPLAY="$DISPLAY" \
    -v "$PROJECT_ROOT:/ws" \
    arg/cg2o
```

---

## Build Library + Examples

```bash
mkdir -p /ws/build && cd /ws/build && \
rm -rf CMakeFiles && rm -f CMakeCache.txt && \
cmake \
  -DCG2O_BUILD_EXAMPLES=ON \
  -DUSE_G2O_SOLVERS=ON \
  -DEXAMPLE_1_OPTIMIZER=ISPD \
  -DEXAMPLE_2_OPTIMIZER=AL \
  -DEXAMPLE_3_OPTIMIZER=BIPM \
  -DCG2O_BUILD_PARDISO=OFF \
  -DCG2O_BUILD_EIGEN_SOLVERS=ON \
  .. && \
make -j$(nproc)
```

```markdown
## To run the examples

After building, you can run the example executables directly from the build directory:

```bash
/ws/build/cg2o/examples/example_1 150 0 0
/ws/build/cg2o/examples/example_2
/ws/build/cg2o/examples/example_3

---

## Standalone MPC Project

```markdown
### Build the MPC example with ISPD (default recommended mode)

This is the main mode used for the proposed solver in the paper. It uses:

- infeasible-start initialization
- the `ISPD` optimizer

```bash
cd /ws/app/cpp_code && \
mkdir -p build && cd build && \
rm -rf CMakeFiles && rm -f CMakeCache.txt && \
cmake .. \
  -DMPC_FEASIBLE_INITIALIZATION=OFF \
  -DMPC_CG2O_OPTIMIZER=ISPD \
  -DMPC_USE_NUMERICAL_JACOBIANS=OFF \
  -DCG2O_DEBUG_LINEAR_SOLVER=OFF && \
make -j$(nproc)

### Build the MPC example with BIPM

If you want to run the MPC application with the **feasible-start Barrier Interior Point Method (BIPM)** instead of the default **ISPD** solver, configure the build as follows.

This mode requires: feasible initialization of the MPC problem
 

Use this configuration when you want to reproduce the feasible-start barrier formulation instead of the infeasible-start primal-dual method.

```bash
cd /ws/app/cpp_code && \
mkdir -p build && cd build && \
rm -rf CMakeFiles && rm -f CMakeCache.txt && \
cmake .. \
  -DMPC_FEASIBLE_INITIALIZATION=ON \
  -DMPC_CG2O_OPTIMIZER=BIPM \
  -DMPC_USE_NUMERICAL_JACOBIANS=OFF \
  -DCG2O_DEBUG_LINEAR_SOLVER=OFF && \
make -j$(nproc)


### Run MPC

```bash
/ws/app/cpp_code/results/data/mpc_g2o_exc
```

### Config File

```bash
nano /ws/app/cpp_code/config/config.yaml
```

 
---

## ROS2 Package

### Build

```bash
source /opt/ros/jazzy/setup.bash
cd /ws/app/ros_ws
rm -rf build install log
colcon build
```

### Run To Test

```bash
        source  /ws/app/ros_ws/install/setup.bash
        
	ros2 run acc_control_mpc_g2o acc_control_mpc_g2o_ros --ros-args -p numberOfIterations:=40 -p horizon_length:=3 & \
	sleep 2; \
	ros2 action send_goal /acc_control_mpc_g2o_action acc_interfaces/action/AccControlMPC "{horizon_length: 3, v_p: 15.0, v_h: 12.0, a_p: 0.5, a_p_weighted: 0.5, d_h: 14.5, force_prev: 10.0}";\
	pkill -9 -f acc_control_mpc_g2o_ros
```

---

## Sparse Linear Solvers

Supports:

- CHOLMOD
- CSparse
- Eigen Sparse Solver*
- UMFPACK LU*
- PARDISO LU*

* Useful for constrained KKT systems.

---

## Debugging

Enable:

```bash
-DCG2O_DEBUG_LINEAR_SOLVER=ON
```

Exports linear systems to text files.

---

## Minimal Example

```cpp
cg2o::SparseOptimizerISPD optimizer;

optimizer.addVertex(...);
optimizer.addEdge(...);
optimizer.addEdgeEq(...);
optimizer.addEdgeIneq(...);

optimizer.initializeOptimization();
optimizer.optimize(50);
```

---

# 
---

## Citation

```bibtex
@article{cg2o2026,
  title={cg2o: A Constrained Factor Graph Solver via the Infeasible Start Interior Point Method for Nonlinear MPC},
  journal={},
  year={ }
}
```

---

 
---

## Final Note

If your problem has:

- sparse structure
- nonlinear constraints
- real-time requirements
- graph structure

then **cg2o** was built for it.
