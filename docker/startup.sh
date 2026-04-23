#!/bin/bash

# 1. Setup the shell environment
if [ -f /ws/scripts/bashSetup.sh ]; then
    source /ws/scripts/bashSetup.sh
fi

# 2. Display the instructions
if [ -f /ws/scripts/welcomeMessage.sh ]; then
    bash /ws/scripts/welcomeMessage.sh
fi

# 3. Stay in the shell
cd /ws
exec /bin/bash
