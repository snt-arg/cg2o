cmake_minimum_required(VERSION 3.5)
project(acc_config)

find_package(ament_cmake REQUIRED)



install(
  DIRECTORY acc_controller/
  DESTINATION share/${PROJECT_NAME}/acc_controller/
)


ament_package()
