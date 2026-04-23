# 🚗 ACC MPC Controller (g2o-based)

This package implements an Adaptive Cruise Control (ACC) Model Predictive Controller using the `cg2o` optimization library.  
The controller is exposed as a ROS 2 **action server**.

---

## 📦 Requirements

- ROS 2 **Jazzy**
- `colcon` build tool
- Dependencies installed (Eigen, g2o, cg2o, yaml-cpp, etc.)

---

## ⚙️ Build Instructions

From your workspace (`ros_ws`):

```bash
source /opt/ros/jazzy/setup.bash
colcon build
source install/setup.bash
```


---

## ▶️ Run the Node

```bash
ros2 run acc_control_mpc_g2o acc_control_mpc_g2o_ros \
  --ros-args \
  -p numberOfIterations:=120 \
  -p algIneqSettings:="[0.01, 1500.0, 10.0]"
```

### Parameters

- `numberOfIterations` (int): number of optimizer iterations  
- `algIneqSettings` (double array):  
  `[kappa_initial, kappa_final, update_factor]`

> ⚠️ All values in `algIneqSettings` must be floating-point numbers.

---

## 📡 Send an Action Request

Open every new terminal for the workspace (ros_ws):

```bash
source /opt/ros/jazzy/setup.bash
source install/setup.bash
```
run the service:

 
```bash
ros2 run acc_control_mpc_g2o acc_control_mpc_g2o_ros --ros-args -p numberOfIterations:=120 -p horizon_length:=4

```



Send a goal from another terminal:
```bash
ros2 action send_goal \
  /acc_control_mpc_g2o_action \
  acc_interfaces/action/AccControlMPC \
  "{horizon_length: 4, v_p: 0.0, v_h: 0.216081, a_p: 0.0, a_p_weighted: 0.5, d_h: 10.0, force_prev: 1850.0}"
```

Note: The horizon length in the sent goal should match the N configuration value in "/ws/app/ros_ws/packages/acc_stuff/acc_config/acc_controller" if the argument -p horizon_length:=x is not given when run the nodeusing: ros2 run acc_control_mpc_g2o acc_control_mpc_g2o_ros.

---

## 🧠 What Happens

- The node runs an MPC controller using `cg2o`
- It receives goals via ROS 2 actions
- It computes optimal control inputs
- It returns the result as an action response

---
