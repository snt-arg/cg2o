# Build and Run

Before building this project, make sure both `g2o` and `cg2o` are installed.

If `cg2o` is not installed, follow the installation instructions in the README file located in the root directory of this project.

---

## Build

Create the build directory:

```bash
mkdir build
cd build
```

# Option 1: ISPD Optimizer
```bash
cmake .. \
  -DMPC_FEASIBLE_INITIALIZATION=OFF \
  -DMPC_CG2O_OPTIMIZER=ISPD \
  -DMPC_USE_NUMERICAL_JACOBIANS=OFF
```


# Option 2: BIPM Optimizer
```bash
cmake .. \
  -DMPC_FEASIBLE_INITIALIZATION=ON \
  -DMPC_CG2O_OPTIMIZER=BIPM \
  -DMPC_USE_NUMERICAL_JACOBIANS=OFF
```


make -j

## Run
# From the build directory:
```bash
taskset -c 0-11 ../results/data/mpc_g2o_exc
```

 
