# cg2o: Constrained Factor Graph Optimization

> **cg2o** (Constrained General Graph Optimization) extends **g2o** with nonlinear **equality and inequality
> constraints** as native factors, solved by an **infeasible-start primal–dual interior-point method (ISPD-IPM)**
> that keeps the sparse structure of the factor graph. It targets nonlinear MPC, trajectory optimization and other
> constrained sparse problems in robotics.

This repository contains the library of the paper

**cg2o: Constrained Factor Graph Optimization via an Infeasible-Start Primal–Dual Interior-Point Method**

---

## Important notes

* **MACC benchmark (paper, Sec. VI).** The code and scripts that reproduce the adaptive cruise control results of
  the paper (Tables IV, V, VII and Figs. 2, 3, 4, 7) are in a separate repository:
  **[MACC benchmark of cg2o](https://anonymous.4open.science/r/cg2o-macc-benchmark/README.md)**. It builds this library at the commit used for the paper and
  writes the paper's tables and figures from every run.
* **Experiment video.** The video of the paper is in the MACC benchmark repository: **`paper/Video.mp4`**.
  It is also available on Google Drive:
  [Experiment Video](https://drive.google.com/file/d/1D5Szbq0rx6m6J8egvvjbGFpykSk1dnAy/view)
  (the link printed in the paper may not open; use this one).
* **Quadrotor benchmark (paper, Sec. VII).** the code of the quadrotor MPC benchmark (PX4
  software-in-the-loop flights) will be linked here: **[QUADROTOR_BENCHMARK_URL](https://anonymous.4open.science/r/cg2o-quadrotor-benchmark-FCC4/README.md)**. Its videos will be placed in
  [`videos/`](videos/).

---

## Overview

Factor graph optimization is a standard tool in robotics for SLAM, state estimation, calibration and trajectory
optimization. Classical frameworks such as g2o solve unconstrained nonlinear least-squares problems; constraints
have to be added as soft penalties.

**cg2o** keeps the factor-graph modelling of g2o and adds hard constraints:

* equality constraints `h(x) = 0` (e.g. system dynamics, kinematic or calibration constraints);
* inequality constraints `g(x) <= 0` (e.g. actuator limits, speed bounds, safety distances, collision avoidance).

Each constraint is a factor (edge) of the graph, so existing g2o vertices, cost edges and models can be reused.

Results reported in the paper:

* adaptive cruise control: up to **54 % fewer iterations** than a tuned barrier-based factor graph solver from the
  same feasible start; from an infeasible cold start **4.6 to 11.2 times faster** than tuned IPOPT;
* quadrotor MPC: all 20 PX4 software-in-the-loop flights completed from randomized starts, a feasible solution in
  every solve despite 91 % obstacle-infeasible initializations, **3.8×** (cold) and **6.7×** (warm) faster than
  IPOPT with HSL solvers.

---

## Features

**Solvers**

| class | method | start |
|---|---|---|
| `cg2o::SparseOptimizerISPD` | **ISPD-IPM**, infeasible-start primal–dual interior-point method (main method of the paper, recommended) | any point, also infeasible |
| `cg2o::SparseOptimizerBIPM` | barrier interior-point method | strictly feasible point |
| `cg2o::SparseOptimizerAL` | augmented Lagrangian method | any point |

**ISPD-IPM**

* infeasible initialization: no feasible starting point needed;
* a single optimization loop instead of nested outer / inner loops;
* inequality constraints as dedicated factors derived from the KKT conditions, with auxiliary (slack) variables;
* lower-bounded adaptive barrier update;
* minimum-threshold backtracking with positivity-preserving safeguards for the auxiliary variables;
* receding-horizon warm starts of primal and dual variables for MPC;
* all options available through setters (`setKappaInitial`, `setKappaFinal`, `setKappaUpdateFactor`, `setTau`,
  `setIneqBacktrackingStepMin`, `setAuxBacktrackingStepMin`, `setAuxScalingFactor`, `setEqualityRegularization`, ...).

**Linear solvers for the indefinite KKT system**

| family | solvers | notes |
|---|---|---|
| unregularized LDLᵀ | **QDLDLT** (`LinearSolverEigenQDLDLT`) | always built; exact KKT step with iterative refinement and LU fallback |
| sparse LU | UMFPACK (`LinearSolverEigenUmfPackLU`), Eigen LU, PARDISO (`LinearSolverEigenPardisoLU`, needs Intel MKL) | UMFPACK is built by default |
| regularized LDLᵀ / Cholesky | CHOLMOD, Eigen LDLT, CHOLMOD fast | need a small equality regularization (`setEqualityRegularization`, e.g. 1e-8) |
| others (Eigen) | QR, SVD, dense LU / QR, BiCGSTAB, LSCG | for testing and small problems |

**Other**

* works with the g2o vertex / edge / block-solver API;
* installs a CMake package (`find_package(cg2o)`, targets `cg2o::...`);
* optional export of every linear system for debugging (`-DCG2O_DEBUG_LINEAR_SOLVER=ON`);
* Docker image with all dependencies, and a ROS 2 example.

---

## Installation

Tested on Ubuntu 22.04 and 24.04 (GCC, CMake ≥ 3.22, C++20).

### 1. Dependencies

```bash
sudo apt update
sudo apt install -y git build-essential cmake libeigen3-dev libfmt-dev libsuitesparse-dev libomp-dev
```

### 2. Build and install

```bash
git clone https://github.com/snt-arg/cg2o.git
cd cg2o
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DUSE_G2O_SOLVERS=ON -DCG2O_BUILD_EIGEN_SOLVERS=ON \
      -DCMAKE_INSTALL_PREFIX="$PWD/install/cg2o"
cmake --build build -j"$(nproc)"
cmake --install build
```

If g2o is not installed, the configure step downloads and builds the tested g2o revision itself
(`scripts/build_g2o.sh`, into `install/externals/g2o/`). An installed g2o is used instead if CMake finds it; it must
be built with `G2O_NO_IMPLICIT_OWNERSHIP_OF_OBJECTS=ON` and with the same `-march=native` setting as cg2o.

### 3. Build options

| option | default | meaning |
|---|---|---|
| `USE_G2O_SOLVERS` | OFF | g2o's CHOLMOD and CSparse linear solvers |
| `CG2O_BUILD_UMFPACK` | ON | UMFPACK LU solver |
| `CG2O_BUILD_EIGEN_SOLVERS` | OFF | Eigen solvers (LDLT, LU, QR, ...) |
| `CG2O_BUILD_PARDISO` | OFF | PARDISO LU solver (needs Intel oneAPI MKL in `/opt/intel/oneapi/mkl/latest`) |
| `CG2O_BUILD_EXAMPLES` | OFF | the examples in `cg2o/examples/` |
| `CG2O_MARCH_NATIVE` | ON | `-march=native`; every program that links cg2o must use the same setting |
| `CG2O_DEBUG_LINEAR_SOLVER` | OFF | write every linear system to text files |

### 4. Use cg2o in your project

```cmake
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")   # same setting as cg2o (CG2O_MARCH_NATIVE)
find_package(cg2o REQUIRED)
target_link_libraries(my_app PRIVATE cg2o::cg2o)
```

Point CMake to cg2o and, if cg2o built it, to its g2o:

```bash
cmake -S . -B build \
  -DCMAKE_PREFIX_PATH="/path/to/cg2o/install/cg2o;/path/to/cg2o/install/externals/g2o/Release"
```

### 5. Examples

```bash
cmake -S . -B build -DCG2O_BUILD_EXAMPLES=ON -DUSE_G2O_SOLVERS=ON -DCG2O_BUILD_EIGEN_SOLVERS=ON
cmake --build build -j"$(nproc)"
./build/cg2o/examples/example_1 150 0 0
./build/cg2o/examples/example_2
./build/cg2o/examples/example_3
```

The optimizer of each example can be chosen at configure time, e.g. `-DEXAMPLE_1_OPTIMIZER=ISPD`,
`-DEXAMPLE_2_OPTIMIZER=AL`, `-DEXAMPLE_3_OPTIMIZER=BIPM`.

---

## Minimal example

```cpp
cg2o::SparseOptimizerISPD optimizer;     // or SparseOptimizerBIPM / SparseOptimizerAL

optimizer.addVertex(v);                  // g2o vertices
optimizer.addEdge(cost_edge);            // g2o cost edges
optimizer.addEdgeEq(eq_edge);            // cg2o::BaseFixedSizedEdgeEq   : h(x) = 0
optimizer.addEdgeIneq(ineq_edge);        // cg2o::BaseFixedSizedEdgeIneq : g(x) <= 0

optimizer.initializeOptimization();
optimizer.optimize(50);
```

---

## Build with Docker

### Build the image

```bash
docker build \
  --build-arg ENABLE_PARDISO=OFF \
  -t arg/cg2o \
  -f ./docker/Dockerfile .
```

### Run the container

```bash
./scripts/runDocker.sh
```

The script mounts the workspace, forwards the display and sets the user permissions. Manually:

```bash
docker run -it --rm \
  --name cg2o \
  --user $(id -u):$(id -g) \
  -v /tmp/.X11-unix:/tmp/.X11-unix \
  -e DISPLAY="$DISPLAY" \
  -v "$(pwd):/ws" \
  arg/cg2o
```

Run the commands from the project root. A second shell in the running container: `docker exec -it cg2o bash`.

Inside the container, build the library and the examples:

```bash
mkdir -p /ws/build && cd /ws/build && \
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

---

## Adaptive cruise control application

`app/adaptive_cruise_control/cpp_code` is a standalone MPC application that uses the installed library. (The
reproduction of the paper's MACC results is the separate [MACC benchmark](MACC_BENCHMARK_URL) repository.)

With ISPD (recommended):

```bash
cd /ws/app/adaptive_cruise_control/cpp_code && \
mkdir -p build && cd build && \
cmake \
  -DMPC_FEASIBLE_INITIALIZATION=OFF \
  -DMPC_CG2O_OPTIMIZER=ISPD \
  -DMPC_USE_NUMERICAL_JACOBIANS=OFF \
  -DCG2O_DEBUG_LINEAR_SOLVER=OFF \
  .. && \
make -j$(nproc)
/ws/app/adaptive_cruise_control/cpp_code/results/data/mpc_g2o_exc
```

With the barrier method (BIPM), which needs a feasible initialization:

```bash
cmake .. \
  -DMPC_FEASIBLE_INITIALIZATION=ON \
  -DMPC_CG2O_OPTIMIZER=BIPM \
  -DMPC_USE_NUMERICAL_JACOBIANS=OFF \
  -DCG2O_DEBUG_LINEAR_SOLVER=OFF && \
make -j$(nproc)
```

Controller, solver and scenario settings: `app/adaptive_cruise_control/cpp_code/config/config.yaml`.

---

## ROS 2 package

`app/adaptive_cruise_control/ros_ws` contains an adaptive cruise control MPC node implemented with cg2o, as a
reference for integrating cg2o into ROS 2.

```bash
source /opt/ros/jazzy/setup.bash && \
cd /ws/app/adaptive_cruise_control/ros_ws && \
rm -rf build install log && \
colcon build
```

Run the node and send a test goal:

```bash
bash -lc '
source /ws/app/adaptive_cruise_control/ros_ws/install/setup.bash

ros2 run acc_control_mpc_g2o acc_control_mpc_g2o_ros \
    --ros-args \
    -p numberOfIterations:=40 \
    -p horizon_length:=3 &

until ros2 action list -t | grep -q "/acc_control_mpc_g2o_action.*acc_interfaces/action/AccControlMPC"; do
    sleep 0.5
done

ros2 action send_goal \
    /acc_control_mpc_g2o_action \
    acc_interfaces/action/AccControlMPC \
    "{horizon_length: 3, v_p: 15.0, v_h: 12.0, a_p: 0.5, a_p_weighted: 0.5, d_h: 14.5, force_prev: 10.0}"

pkill -INT -f acc_control_mpc_g2o_ros
'
```

`numberOfIterations` is the iteration limit of the optimizer, `horizon_length` the MPC horizon.

---

## Citation

```bibtex
@article{cg2o2026,
  title   = {cg2o: Constrained Factor Graph Optimization via an Infeasible-Start Primal--Dual Interior-Point Method},
  journal = {},
  year    = {2026}
}
```
