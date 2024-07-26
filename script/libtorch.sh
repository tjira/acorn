#!/bin/bash

# download and extract libtorch
wget https://github.com/pytorch/pytorch/releases/download/v2.4.0/pytorch-v2.4.0.tar.gz && tar -xzvf pytorch-v2.4.0.tar.gz && mv pytorch-v2.4.0 libtorch && rm pytorch-v2.4.0.tar.gz

# configure libtorch
cd libtorch && cmake -B build -DBUILD_PYTHON=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install" -DUSE_FBGEMM=OFF && cd ..

# build and install libtorch
cd libtorch && cmake --build build --parallel 2 && cmake --install build && cd ..

# pip install pyyaml typing-extensions
