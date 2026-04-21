#!/bin/bash

# Either  bashSetup.sh --clean or bashSetup.sh
BASHRC="$HOME/.bashrc"
MODE=$1  # Takes the first argument (e.g., --clean)

# --- Helper Functions ---

add_to_bashrc() {
    grep -qF "$1" "$BASHRC" || echo "$1" >> "$BASHRC"
}

remove_from_bashrc() {
    # Escape the string for sed and delete matching lines
    sed -i "\|$1|d" "$BASHRC"
}

# --- Management Logic ---

process_line() {
    local LINE=$1
    if [ "$MODE" == "--clean" ]; then
        remove_from_bashrc "$LINE"
    else
        add_to_bashrc "$LINE"
    fi
}

#------------------------------ What to add
# Add ROS and Workspace sourcing
process_line "source /opt/ros/jazzy/setup.bash"
process_line "if [ -f /ws/app/ros_ws/install/setup.bash ]; then source /ws/app/ros_ws/install/setup.bash; fi"

# Add Aliases
process_line 'alias build_all="python3 /ws/scripts/build_all.py"'
process_line 'alias run_all="python3 /ws/scripts/run_all.sh"'

# Add Threads
#process_line "export OPENBLAS_NUM_THREADS=8"
#process_line "export OMP_NUM_THREADS=8"
#process_line "export MKL_NUM_THREADS=8"

source "$BASHRC"

# --- Finalize ---

if [ "$MODE" == "--clean" ]; then
    echo "Cleanup complete. Your .bashrc is now clean of CG2O entries."
else
    echo "Setup applied. Aliases and ROS environment are ready."
    source "$BASHRC"
fi
