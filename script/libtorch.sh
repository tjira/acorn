#!/bin/bash

# download libtorch
mkdir -p external && wget -O external/libtorch.tar.gz https://github.com/pytorch/pytorch/releases/download/v2.4.0/pytorch-v2.4.0.tar.gz

# unpack libtorch
cd external && rm -rf libtorch && tar -xzvf libtorch* && mv pytorch* libtorch && cd -

# configure libtorch
cd external/libtorch && cmake -B build -DBUILD_PYTHON=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install" -DUSE_DISTRIBUTED=OFF -DUSE_FBGEMM=OFF -DUSE_OPENMP=OFF && cd -

# compile and install libtorch
cd external/libtorch && cmake --build build --parallel 2 && cmake --install build && cd -

# copy the compiled library
cp -r external/libtorch/install/include external/libtorch/install/lib external/

# remove redundant files
rm -rf external/libtorch external/lib/cmake external/lib/pkgconfig external/lib/libtorch_global_deps.so external/lib/*.a
