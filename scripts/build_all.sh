#!/bin/bash

set -e   # stop on error

echo "🔧 Building ROS workspace..."

source /opt/ros/jazzy/setup.bash

cd /ws/app/ros_ws

rm -rf build install log

colcon build

echo "✅ ROS workspace built"
echo

echo "🔧 Building cpp_code..."

cd /ws/app/cpp_code

rm -rf build
mkdir -p build
cd build

cmake ..
make -j"$(nproc)"

echo "✅ cpp_code built"


source ~/.bashrc
