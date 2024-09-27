#!/bin/bash

# download libeigen
mkdir -p external && git clone https://gitlab.com/libeigen/eigen.git external/libeigen

# configure libeigen
cd external/libeigen && cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX="$PWD/install" -DEIGEN_BUILD_BLAS=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_LAPACK=OFF -DEIGEN_BUILD_TESTING=OFF -DEIGEN_TEST_NOQT=ON && cd -

# compile and install libeigen
cd external/libeigen && cmake --build build --parallel 2 && cmake --install build && cd -

# copy the compiled library
cp -r external/libeigen/install/include external/
