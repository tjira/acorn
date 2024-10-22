#!/bin/bash

CORES=2

# set paths
export CPLUS_INCLUDE_PATH="$PWD/external/include:$CPLUS_INCLUDE_PATH"; export LD_LIBRARY_PATH="$PWD/external/lib:$LD_LIBRARY_PATH"; export LIBRARY_PATH="$PWD/external/lib:$LIBRARY_PATH"

# make the folders
mkdir -p external && mkdir -p external/include mkdir && -p external/lib

# download argparse
wget -O external/include/argparse.hpp https://raw.githubusercontent.com/p-ranav/argparse/master/include/argparse/argparse.hpp

# download json
wget -O external/include/json.hpp https://raw.githubusercontent.com/nlohmann/json/develop/single_include/nlohmann/json.hpp

# download exprtk
wget -O external/include/exprtk.hpp https://raw.githubusercontent.com/ArashPartow/exprtk/master/exprtk.hpp

# download fftw
wget -O external/libfftw.tar.gz https://www.fftw.org/fftw-3.3.10.tar.gz

# download libint
wget -O external/libint.tar.gz https://github.com/evaleev/libint/releases/download/v2.9.0/libint-2.9.0-mpqc4.tgz

# download libnuma
wget -O external/libnuma.tar.gz https://github.com/numactl/numactl/releases/download/v2.0.18/numactl-2.0.18.tar.gz

# download libopenblas
wget -O external/openblas.tar.gz https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.28/OpenBLAS-0.3.28.tar.gz

# download libtorch
wget -O external/libtorch.tar.gz https://github.com/pytorch/pytorch/releases/download/v2.5.0/pytorch-v2.5.0.tar.gz

# unpack the archives
cd external && for ARCHIVE in *.tar.gz; do tar -xzf $ARCHIVE --warning=no-unknown-keyword; done; cd ..

# compile fftw
cd external/fftw-3.3.10 && cmake -B build \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
&& cmake --build build --parallel $CORES && cmake --install build && cp -r install/* .. ; cd ../..

# compile libint
cd external/libint-2.9.0 && cmake -B build \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
&& cmake --build build --parallel $CORES && cmake --install build && cp -r install/* .. ; cd ../..

# compile numa
cd external/numactl-2.0.18 && ./configure \
    --enable-shared="no" \
    --enable-static="yes" \
    --prefix="$PWD/install" \
&& make -j$CORES && make install && cp -r install/* .. ; cd ../..

# compile openblas
cd external/OpenBLAS-0.3.28 && cmake -B build \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_FLAGS="-Wno-error=incompatible-pointer-types" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
    -DNOFORTRAN=ON \
&& cmake --build build --parallel $CORES && cmake --install build && cp -r install/* .. ; cd ../..

# compile libtorch
cd external/pytorch-v2.5.0 && cmake -B build \
    -DBUILD_PYTHON=OFF \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
    -DUSE_CUDA=OFF \
    -DUSE_DISTRIBUTED=OFF \
    -DUSE_FBGEMM=OFF \
    -DUSE_MKLDNN=OFF \
    -DUSE_OPENMP=OFF \
    -DUSE_PYTORCH_QNNPACK=OFF \
    -DUSE_ROCM=OFF \
    -DUSE_XNNPACK=OFF \
    -DUSE_XPU=OFF \
&& cmake --build build --parallel $CORES && cmake --install build && cp -r build/lib .. && cp -r install/* .. ; cd ../..
