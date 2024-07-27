#!/bin/bash

# download libtorch
git clone --depth 1 --recursive https://github.com/pytorch/pytorch.git libtorch

# configure libtorch
cd libtorch && cmake -B build -DBUILD_PYTHON=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install" -DUSE_DISTRIBUTED=OFF -DUSE_FBGEMM=OFF -DUSE_OPENMP=OFF && cd ..

# build and install libtorch
cd libtorch && cmake --build build --parallel 2 && cmake --install build && cd ..
