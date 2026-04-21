FROM docker.io/ros:jazzy

ENV DEBIAN_FRONTEND=noninteractive
ARG G2O_COMMIT=eec325a1
ARG CG2O_COMMIT=eec325a1

ARG ENABLE_PARDISO=OFF
ENV ENABLE_PARDISO=${ENABLE_PARDISO}


# -------------------------------------------------
# System dependencies
# -------------------------------------------------
RUN apt-get update && apt-get install -y \
    nano \
    vim \
    trash-cli \
    git \
    wget \
    cmake \
    build-essential \
    clang \
    gdb \
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    python3-matplotlib \
    libeigen3-dev \
    libspdlog-dev \
    libsuitesparse-dev \
    libmetis-dev \
    libfmt-dev \
    libomp-dev \
    qtdeclarative5-dev \
    qt5-qmake \
    libqglviewer-dev-qt5 \
    python3-colcon-common-extensions \
    python3-rosdep \
    ros-jazzy-common-interfaces \
    ros-jazzy-rviz2 \
    ros-jazzy-pcl-ros \
    ros-jazzy-cv-bridge \
    ros-jazzy-image-transport-plugins \
    && rm -rf /var/lib/apt/lists/*

# Set clang as default compiler
RUN update-alternatives --install /usr/bin/c++ c++ /usr/bin/clang++ 100 && \
    update-alternatives --install /usr/bin/cc cc /usr/bin/clang 100

# -------------------------------------------------
# Build and install g2o
# -------------------------------------------------
RUN cd /tmp && \
    git clone https://github.com/RainerKuemmerle/g2o.git && \
    cd g2o && \
    git checkout ${G2O_COMMIT} && \
    sed -i 's/^  void update(const double\* update);/  virtual void update(const double* update);/' \
        g2o/core/sparse_optimizer.h && \
    mkdir build && cd build && \
    cmake \
      -DCMAKE_BUILD_TYPE=Release \
      -DG2O_USE_OPENMP=ON \
      -DG2O_BUILD_APPS=OFF \
      -DG2O_BUILD_EXAMPLES=OFF \
      -DG2O_USE_CHOLMOD=ON \
      -DG2O_USE_CSPARSE=ON \
      -DCMAKE_CXX_FLAGS="-DEIGEN_USE_THREADS -DEIGEN_USE_OPENMP -fopenmp" \
      .. && \
    make -j"$(nproc)" && \
    make install && \
    ldconfig && \
    rm -rf /tmp/g2o

# -------------------------------------------------
# Optional Intel oneAPI MKL for PARDISO
# -------------------------------------------------
RUN if [ "$ENABLE_PARDISO" = "ON" ]; then \
      cd /tmp && \
      wget -q https://registrationcenter-download.intel.com/akdlm/IRC_NAS/d0df6732-bf5c-493b-a484-6094bea53787/intel-oneapi-hpc-toolkit-2025.1.0.666_offline.sh && \
      chmod +x intel-oneapi-hpc-toolkit-2025.1.0.666_offline.sh && \
      sh ./intel-oneapi-hpc-toolkit-2025.1.0.666_offline.sh -a --silent --eula accept && \
      rm -f intel-oneapi-hpc-toolkit-2025.1.0.666_offline.sh; \
    else \
      echo "Skipping Intel oneAPI installation"; \
    fi

# -------------------------------------------------
# Environment
# -------------------------------------------------
ENV LD_LIBRARY_PATH=/usr/local/lib:/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH}
ENV PYTHONPATH=/usr/local/lib/python3.12/dist-packages:${PYTHONPATH}

# -------------------------------------------------
# Build and install cg2o
# -------------------------------------------------
RUN cd /tmp && \
    git clone https://github.com/snt-arg/cg2o.git && \
    cd cg2o && \
    mkdir build && cd build && \
    cmake \
       -DCG2O_BUILD_EXAMPLES=OFF \
       -DUSE_G2O_SOLVERS=ON \
       -DCG2O_BUILD_PARDISO=${ENABLE_PARDISO} \
       -DCG2O_BUILD_EIGEN_SOLVERS=ON \
       .. && \
    make -j"$(nproc)" && \
    make install && \
    ldconfig && \
    rm -rf /tmp/cg2o


# -------------------------------------------------
# User setup
# -------------------------------------------------
RUN useradd -s /bin/bash -m cg2o && \
    echo "cg2o:cg2o" | chpasswd && \
    adduser cg2o sudo && \
    mkdir -p /etc/sudoers.d && \
    echo "ALL ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/anyuser && \
    chmod 440 /etc/sudoers.d/anyuser

RUN mkdir -p /ws && chmod 777 /ws

# -------------------------------------------------
# Startup
# -------------------------------------------------
COPY docker/startup.sh /usr/local/bin/startup.sh
RUN chmod +x /usr/local/bin/startup.sh

USER cg2o
WORKDIR /ws
ENTRYPOINT ["/usr/local/bin/startup.sh"]
