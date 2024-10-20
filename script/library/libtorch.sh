#!/bin/bash

# download the library
mkdir -p external && wget -O external/libtorch.tar.gz https://github.com/pytorch/pytorch/releases/download/v2.5.0/pytorch-v2.5.0.tar.gz && cd external && rm -rf libtorch && tar -xzvf libtorch* && mv pytorch* libtorch && cd -

# configure the library
cd external/libtorch && cmake -B build -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install" -DUSE_CUDA=OFF && cd -

# compile and install the library
cd external/libtorch && cmake --build build --parallel 2 && cmake --install build && cd -

# copy the compiled library
cp -r external/libtorch/install/include external/libtorch/install/lib external/
