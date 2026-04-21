#!/bin/bash

# 1. Setup the shell environment
if [ -f /ws/docker/bashSetup.sh ]; then
    source /ws/docker/bashSetup.sh
fi

# 2. Display the instructions
if [ -f /ws/docker/welcomeMessage.sh ]; then
    bash /ws/docker/welcomeMessage.sh
fi

# 3. Stay in the shell
cd /ws
exec /bin/bash
