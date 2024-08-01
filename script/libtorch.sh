#!/bin/bash

# download libtorch
mkdir -p external && wget -O external/libtorch.tar.gz https://github.com/pytorch/pytorch/releases/download/v2.4.0/pytorch-v2.4.0.tar.gz

# unpack libtorch
cd external && rm -rf libtorch && tar -xzvf libtorch* && mv pytorch* libtorch && cd -

# configure libtorch
cd external/libtorch && cmake -B build -DBUILD_PYTHON=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install" -DUSE_CUDA=OFF -DUSE_DISTRIBUTED=OFF -DUSE_FBGEMM=OFF -DUSE_MAGMA=OFF -DUSE_MKLDNN=OFF -DUSE_NNPACK=OFF -DUSE_NUMPY=OFF -DUSE_OPENMP=OFF -DUSE_PYTORCH_QNNPACK=OFF -DUSE_ROCM=OFF -DUSE_XNNPACK=OFF -DUSE_XPU=OFF && cd -

# compile and install libtorch
cd external/libtorch && cmake --build build --parallel 2 && cmake --install build && cd -

# remove static libraries and copy the compiled library
rm external/libtorch/install/lib/*.a && cp -r external/libtorch/install/include external/libtorch/install/lib external/

# remove redundant files
rm -rf external/libtorch external/lib/cmake external/lib/pkgconfig external/lib/libtorch_global_deps.so
