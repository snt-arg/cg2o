#!/bin/bash

  docker build --build-arg ENABLE_PARDISO=OFF -t arg/cg2o -f ../docker/Dockerfile ..
