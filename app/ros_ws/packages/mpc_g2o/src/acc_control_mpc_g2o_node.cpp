#include "acc_interfaces/srv/empty.hpp"
#include "acc_mpc_control_g2o.h"
#include "acc_prediction.h"
#include <acc_interfaces/action/acc_control_mpc.hpp>
#include <ament_index_cpp/get_package_share_directory.hpp>
#include <chrono>
#include <filesystem> // Added for filesystem operations
#include <memory>
#include <rclcpp/rclcpp.hpp>
#include <rclcpp_action/rclcpp_action.hpp>
#include <sstream>
#include <string>
#include <thread>
#include <yaml-cpp/yaml.h>

class ACCMPCControlG2ONode : public rclcpp::Node {
public:
  using AccControlMPC = acc_interfaces::action::AccControlMPC;
  using Empty = acc_interfaces::srv::Empty;
  using GoalHandleAccControlMPC =
      rclcpp_action::ServerGoalHandle<AccControlMPC>;

  ACCMPCControlG2ONode() : Node("acc_mpc_control_g2o") {
// Initialize log file in workspace/data/logs/acc_control_mpc_g2o/
#ifndef NDEBUG
    RCLCPP_INFO(this->get_logger(), "Debug mode is ON.");
    initializeLogFile();
#endif
    // Get package paths
    std::string pkg_path =
        ament_index_cpp::get_package_share_directory("acc_config");
    std::string config_mpc_parameters_yaml_file_name =
        pkg_path + "/acc_controller/config_acc_mpc_controller.yaml";
    int N_from_yaml =
        getParameterByName(config_mpc_parameters_yaml_file_name, "N");

    // Declare parameters
    this->declare_parameter("numberOfIterations", 30);
    this->declare_parameter("horizon_length", N_from_yaml);

    this->declare_parameter<std::vector<double>>(
        "algIneqSettings",
        std::vector<double>{.5, 1500,
                            8}); // kappa_initial, kappa_final, update_factor
    // Get parameters
    _numberOfIterations = this->get_parameter("numberOfIterations").as_int();
    _mpcHorizon = this->get_parameter("horizon_length").as_int();

    std::vector<double> algIneqSettings =
        this->get_parameter("algIneqSettings").as_double_array();

    RCLCPP_INFO(this->get_logger(), "Horizon Length (N): %d", _mpcHorizon);
    RCLCPP_INFO(this->get_logger(), "numberOfIterations: %d",
                _numberOfIterations);
    RCLCPP_INFO(this->get_logger(), "algIneqSettings: ");
    for (const auto &setting : algIneqSettings) {
      RCLCPP_INFO(this->get_logger(), "%f", setting);
    }

    // Initialize controller and prediction objects
    _controller = std::make_shared<ACCMPCControlG2O>(_numberOfIterations,
                                                     algIneqSettings);
    _controller->initializeOptimizer();
    _controller->initializeMPCGraph(_mpcHorizon);
    _prediction_mpc = std::make_shared<AccPrediction>();

    // Setup action server
    action_server_ = rclcpp_action::create_server<AccControlMPC>(
        this, "acc_control_mpc_g2o_action",
        std::bind(&ACCMPCControlG2ONode::handle_goal, this,
                  std::placeholders::_1, std::placeholders::_2),
        std::bind(&ACCMPCControlG2ONode::handle_cancel, this,
                  std::placeholders::_1),
        std::bind(&ACCMPCControlG2ONode::handle_accepted, this,
                  std::placeholders::_1));

    RCLCPP_INFO(this->get_logger(),
                "acc_control_mpc_g2o_action service is ready.");
  }

  ~ACCMPCControlG2ONode() {
    if (log_file_.is_open()) {
      auto now = std::chrono::system_clock::now();
      auto time_t = std::chrono::system_clock::to_time_t(now);
      log_file_ << std::endl
                << "=======================================" << std::endl;
      log_file_ << "Log ended: "
                << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S")
                << std::endl;
      log_file_ << "Total requests processed: " << action_count_ << std::endl;
      log_file_.close();
      RCLCPP_INFO(this->get_logger(), "Log file closed: %s",
                  log_filename_.c_str());
    }
  }

private:
  int _numberOfIterations = 100;
  int _mpcHorizon;
  int _linearSolverType = 10;
  double v_min_ = -1, v_max_ = 35, road_grad_ = 0;
  rclcpp_action::Server<AccControlMPC>::SharedPtr action_server_;
  rclcpp::Service<Empty>::SharedPtr service_;

  std::shared_ptr<ACCMPCControlG2O> _controller;
  std::shared_ptr<AccPrediction> _prediction_mpc;

  // Logging members
  std::string log_filename_;
  std::ofstream log_file_;
  int action_count_ = 0;

  void initializeLogFile() {
    // Get the workspace path (assuming we're in the install directory)
    std::string install_dir =
        ament_index_cpp::get_package_share_directory("acc_control_mpc_g2o");
    std::filesystem::path pkg_path(install_dir);
    std::filesystem::path workspace_path =
        pkg_path.parent_path().parent_path().parent_path().parent_path();
    std::string log_dir =
        workspace_path.string() + "/data/logs/acc_control_mpc_g2o";
    RCLCPP_INFO(this->get_logger(), "Log directory will be: %s",
                log_dir.c_str());
    // Create directories if they don't exist
    std::filesystem::create_directories(log_dir);

    // Create filename with timestamp
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::ostringstream filename;
    filename << log_dir << "/acc_mpc_g2o_"
             << std::put_time(std::localtime(&time_t), "%Y%m%d_%H%M%S")
             << ".txt";

    log_filename_ = filename.str();
    log_file_.open(log_filename_, std::ios::app);

    if (log_file_.is_open()) {
      // Write simplified header
      auto now = std::chrono::system_clock::now();
      auto time_t = std::chrono::system_clock::to_time_t(now);
      log_file_ << "ACC MPC Control G2O Action Log" << std::endl;
      log_file_ << "Started: "
                << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S")
                << std::endl;
      log_file_ << "File: " << log_filename_ << std::endl;
      log_file_ << "=======================================" << std::endl;
      log_file_ << std::endl;

      log_file_.flush();

      RCLCPP_INFO(this->get_logger(), "Logging to: %s", log_filename_.c_str());
    } else {
      RCLCPP_ERROR(this->get_logger(), "Failed to open log file: %s",
                   log_filename_.c_str());
    }
  }

  void logActionToFile(const std::string &action_type, double v_p, double a_p,
                       double a_p_weighted, double v_h, double d_h,
                       double force_prev, int horizon_length,
                       double result_force = 0.0, int flag_solved = 0,
                       double cost_value = 0.0, int num_iterations = 0,
                       double solve_time_ms = 0.0,
                       double desired_velocity = 0.0,
                       double desired_acceleration = 0.0) {
    if (!log_file_.is_open())
      return;

    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);

    // Format for easy copy-paste to terminal
    if (action_type == "REQUEST") {
      log_file_ << "Request:" << std::endl;
      log_file_ << "ros2 action send_goal /acc_control_mpc_g2o_action "
                   "acc_interfaces/action/AccControlMPC \\"
                << std::endl;
      log_file_ << "\"{horizon_length: " << horizon_length
                << ", v_p: " << std::fixed << std::setprecision(6) << v_p
                << ", v_h: " << v_h << ", a_p: " << a_p
                << ", a_p_weighted: " << a_p_weighted << ", d_h: " << d_h
                << ", force_prev: " << force_prev << "}\"" << std::endl;
      log_file_ << std::endl;
    } else if (action_type == "RESULT") {
      log_file_ << "Response:" << std::endl;
      log_file_ << "result_force: " << std::fixed << std::setprecision(6)
                << result_force << ", flag_solved: " << flag_solved
                << ", cost_value: " << cost_value
                << ", num_iterations: " << num_iterations
                << ", solve_time_ms: " << solve_time_ms
                << ", desired_velocity: " << desired_velocity
                << ", desired_acceleration: " << desired_acceleration
                << std::endl;
      log_file_ << std::endl;
      log_file_ << "=======================================" << std::endl;
      log_file_ << std::endl;
    }

    log_file_.flush();

    if (action_type == "REQUEST") {
      action_count_++;
    }
  }

  int getParameterByName(const std::string &file_name,
                         const std::string &param_name) {
    YAML::Node config = YAML::LoadFile(file_name);
    if (config[param_name]) {
      return config[param_name].as<int>();
    } else {
      RCLCPP_ERROR(this->get_logger(), "Parameter %s not found in %s",
                   param_name.c_str(), file_name.c_str());
      return -1;
    }
  }

  void
  logPredictionData(const std::vector<std::vector<double>> &prediction_data,
                    rclcpp::Logger logger) {
    if (prediction_data.empty()) {
      RCLCPP_ERROR(logger, "Prediction data is empty.");
      return;
    } else {
      std::ostringstream oss;
      oss << "v_p_pred:";
      for (const auto &val : prediction_data[0]) {
        oss << "  " << val;
      }
      oss << " ***" << "road_grad_pred: " << prediction_data[1][0] << " "
          << "v_min_pred: " << prediction_data[2][0] << " "
          << "v_max_pred: " << prediction_data[3][0];
      RCLCPP_INFO(logger, "%s", oss.str().c_str());
    }
  }

  rclcpp_action::GoalResponse
  handle_goal(const rclcpp_action::GoalUUID &uuid,
              std::shared_ptr<const AccControlMPC::Goal> goal) {
    RCLCPP_INFO(this->get_logger(), "Received action goal request.");

// Log the incoming action request
#ifndef NDEBUG
    logActionToFile("REQUEST", goal->v_p, goal->a_p, goal->a_p_weighted,
                    goal->v_h, goal->d_h, goal->force_prev,
                    goal->horizon_length);
#endif

    return rclcpp_action::GoalResponse::ACCEPT_AND_EXECUTE;
  }

  rclcpp_action::CancelResponse
  handle_cancel(const std::shared_ptr<GoalHandleAccControlMPC> goal_handle) {
    RCLCPP_INFO(this->get_logger(), "Received action cancel request.");
    return rclcpp_action::CancelResponse::ACCEPT;
  }

  void
  handle_accepted(const std::shared_ptr<GoalHandleAccControlMPC> goal_handle) {
    std::thread{std::bind(&ACCMPCControlG2ONode::execute, this, goal_handle)}
        .detach();
  }

  void execute(const std::shared_ptr<GoalHandleAccControlMPC> goal_handle) {
    RCLCPP_INFO(this->get_logger(),
                "Processing action request in MPC g2o controller ...");
    const auto goal = goal_handle->get_goal();
    int N = goal->horizon_length;
    const auto &v_p = goal->v_p;
    double a_p = goal->a_p;
    double a_p_weighted = goal->a_p_weighted;
    double v_h = goal->v_h;
    double d_h = goal->d_h;
    double force_prev = goal->force_prev;

    // Log the received parameters
    RCLCPP_INFO(this->get_logger(),
                "Received parameters: v_p: %f, a_p: %f, a_p_weighted: %f, v_h: "
                "%f, d_h: %f, "
                "force_prev: %f, N: %d",
                v_p, a_p, a_p_weighted, v_h, d_h, force_prev, N);

    // 1. Reset graph and reinitialize if horizon changes
    if (_mpcHorizon != N) {
      RCLCPP_INFO(this->get_logger(),
                  "[ERROR] Horizon length changed from %d to %d. Check the "
                  "config parameters for N=%d, but the server get request for "
                  "N=%d. Check again that both are equal.",
                  _mpcHorizon, N, _mpcHorizon, N);
    }

    // 2. Update MPC parameters
    double f_t_prev = 0.0, f_b_prev = 0.0;
    if (force_prev > 0) {
      f_t_prev = force_prev;
    } else {
      f_b_prev = std::abs(force_prev);
    }
    _controller->updateMPCParameters(v_h, d_h, f_t_prev, f_b_prev);

    // 3. Compute prediction data
    auto prediction_data = _prediction_mpc->compute_mpc_prediction(
        N, v_p, a_p, a_p_weighted, v_min_, v_max_, road_grad_);
    _controller->updateMPCPrediction(prediction_data);

    logPredictionData(prediction_data, this->get_logger());
#if CG2O_DEBUG_LINEAR_SOLVER
    _controller->printAll();
#endif
    // 4. Apply the MPC control
    auto start_time = std::chrono::high_resolution_clock::now();
    auto num_iterations =
        _controller->optimizer()->optimize(_numberOfIterations);
    auto solve_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start_time)
            .count();

    // 5. Prepare and send the result
    auto result = std::make_shared<AccControlMPC::Result>();
    result->acc_input_force = _controller->mpc()->getForceInput(0);
    result->flag_solved = 1;
    result->cost_value = _controller->optimizer()->chi2();
    result->num_iterations = num_iterations;
    result->solve_time = solve_time;
    result->desired_velocity = _controller->mpc()->getVelocityDesired();
    result->desired_acceleration =
        _controller->mpc()->getAccelerationDesired(v_h);

// Log the result
#ifndef NDEBUG
    logActionToFile("RESULT", v_p, a_p, a_p_weighted, v_h, d_h, force_prev, N,
                    result->acc_input_force, result->flag_solved,
                    result->cost_value, result->num_iterations,
                    result->solve_time, result->desired_velocity,
                    result->desired_acceleration);
#endif

    goal_handle->succeed(result);
    RCLCPP_INFO(
        this->get_logger(),
        "The results: force_input: %f, num_iteration: %d, solve_time(ms):%f ",
        result->acc_input_force, result->num_iterations, result->solve_time);
  }
};

int main(int argc, char **argv) {
  rclcpp::init(argc, argv);
  auto node = std::make_shared<ACCMPCControlG2ONode>();
  rclcpp::spin(node);
  rclcpp::shutdown();
  return 0;
}