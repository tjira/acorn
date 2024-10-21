#!/bin/bash

# compile parameters
SHARED="OFF"; CORES=2; if [ $# -ne 0 ] && [ $# -ne 2 ]; then echo "ARGUMENTS: SHARED(ON/OFF) CORES(1,2,...)"; exit 1; fi; if [ $# -eq 2 ]; then SHARED=$1; CORES=$2; fi

# download the library
mkdir -p external && wget -O external/libtorch.tar.gz https://github.com/pytorch/pytorch/releases/download/v2.5.0/pytorch-v2.5.0.tar.gz && cd external && rm -rf libtorch && tar -xzvf libtorch* && mv pytorch* libtorch && cd -

# configure the library
cd external/libtorch && cmake -B build -DBUILD_PYTHON=OFF -DBUILD_SHARED_LIBS=$SHARED -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install" -DUSE_CUDA=OFF -DUSE_DISTRIBUTED=OFF -DUSE_FBGEMM=OFF -DUSE_MKLDNN=OFF -DUSE_OPENMP=OFF -DUSE_PYTORCH_QNNPACK=OFF -DUSE_ROCM=OFF -DUSE_XNNPACK=OFF -DUSE_XPU=OFF && cd -

# compile and install the library
cd external/libtorch && cmake --build build --parallel $CORES && cmake --install build && cd -

# copy the compiled library
cp -r external/libtorch/build/lib external/ && cp -r external/libtorch/install/include external/libtorch/install/lib external/
