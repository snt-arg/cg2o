# Sources of the solver start points, compiled into the core library (see
# cg2o/core/CMakeLists.txt).
set(start_point_dir ${CMAKE_CURRENT_LIST_DIR})
set(start_point_source
  ${start_point_dir}/ispd/ispd_start_point.cpp
  ${start_point_dir}/ispd/ispd_strategies.cpp
  ${start_point_dir}/al/al_start_point.cpp
  ${start_point_dir}/bipm/bipm_start_point.cpp
)
set(start_point_header
  ${start_point_dir}/start_point_common.h
  ${start_point_dir}/ispd/ispd_start_point.h
  ${start_point_dir}/ispd/ispd_start_point.hpp
  ${start_point_dir}/ispd/ispd_strategies.h
  ${start_point_dir}/ispd/ispd_start_math.h
  ${start_point_dir}/al/al_start_point.h
  ${start_point_dir}/bipm/bipm_start_point.h
)
