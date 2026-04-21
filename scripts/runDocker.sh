#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# We add --user $(id -u):$(id -g) to match your host permissions
docker run -it --rm \
    --name cg2o \
    --user $(id -u):$(id -g) \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    -e DISPLAY="$DISPLAY" \
    -v "$PROJECT_ROOT:/ws" \
    arg/cg2o
