#!/bin/bash

# check the number of arguments
[ $# -ne 3 ] && echo "USAGE: $0 BUILD_TYPE COMPILE_MODE CORES" && exit 1

# assign the arguments
CORES=$3; MODE=$2; TYPE=${1,,}; [ $MODE == "STATIC" ] && STATIC=1 || STATIC=0; [ $MODE == "SHARED" ] && SHARED=1 || SHARED=0; [ $SHARED -eq 0 ] && [ $STATIC -eq 0 ] && echo "INVALID MODE" && exit 1

echo ${TYPE^}

exit

# set environment paths
export CPLUS_INCLUDE_PATH="$PWD/external/include:$CPLUS_INCLUDE_PATH"; export LIBRARY_PATH="$PWD/external/lib:$LIBRARY_PATH"; export Eigen3_DIR="$PWD/external/share/eigen3/cmake"

# make the folders
mkdir -p external && mkdir -p external/include && mkdir -p external/lib

# download argparse
wget -q -O external/include/argparse.hpp https://raw.githubusercontent.com/p-ranav/argparse/master/include/argparse/argparse.hpp

# download json
wget -q -O external/include/json.hpp https://raw.githubusercontent.com/nlohmann/json/develop/single_include/nlohmann/json.hpp

# download exprtk
wget -q -O external/include/exprtk.hpp https://raw.githubusercontent.com/ArashPartow/exprtk/master/exprtk.hpp

# download eigen
wget -q -O external/libeigen.tar.gz https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz

# download fftw
wget -q -O external/libfftw.tar.gz https://www.fftw.org/fftw-3.3.10.tar.gz

# download libint
wget -q -O external/libint.tar.gz https://github.com/evaleev/libint/releases/download/v2.9.0/libint-2.9.0-mpqc4.tgz

# download libnuma
wget -q -O external/libnuma.tar.gz https://github.com/numactl/numactl/releases/download/v2.0.18/numactl-2.0.18.tar.gz

# download libopenblas
wget -q -O external/openblas.tar.gz https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.28/OpenBLAS-0.3.28.tar.gz

# download libtorch
wget -q -O external/libtorch.tar.gz https://github.com/pytorch/pytorch/releases/download/v2.5.0/pytorch-v2.5.0.tar.gz

# unpack the archives
cd external && for ARCHIVE in *.tar.gz; do tar -xzf $ARCHIVE --warning=no-unknown-keyword; done; cd ..

# compile eigen
cd external/eigen-3.4.0 && cmake -B build \
    -DCMAKE_BUILD_TYPE=${TYPE^} \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
    -DEIGEN_BUILD_DOC=OFF \
    -DEIGEN_TEST_NOQT=ON \
&& cmake --build build --parallel $CORES && cmake --install build && cp -r install/* .. && cd ../..

# compile fftw
cd external/fftw-3.3.10 && cmake -B build \
    -DBUILD_SHARED_LIBS=$SHARED \
    -DCMAKE_BUILD_TYPE=${TYPE^} \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
&& cmake --build build --parallel $CORES && cmake --install build && cp -r install/* .. && cd ../..

# compile libint
cd external/libint-2.9.0 && cmake -B build \
    -DBUILD_SHARED_LIBS=$SHARED \
    -DCMAKE_BUILD_TYPE=${TYPE^} \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
&& cmake --build build --parallel $CORES && cmake --install build && cp -r install/* .. && cd ../..

# compile numa
cd external/numactl-2.0.18 && ./configure \
    --enable-shared=$([ $SHARED == 1 ] && echo "yes" || echo "no") \
    --enable-static=$([ $STATIC == 1 ] && echo "yes" || echo "no") \
    --prefix="$PWD/install" \
&& make -j$CORES && make install && cp -r install/* .. && cd ../..

# compile shared openblas
cd external/OpenBLAS-0.3.28 && cmake -B build \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_BUILD_TYPE=${TYPE^} \
    -DCMAKE_C_FLAGS="-Wno-error=incompatible-pointer-types" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
    -DNOFORTRAN=ON \
&& cmake --build build --parallel $CORES && cmake --install build && cp -r install/* .. && cd ../..

# compile static openblas
[ $STATIC == 1 ] && cd external/OpenBLAS-0.3.28 && cmake -B build \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_BUILD_TYPE=${TYPE^} \
    -DCMAKE_C_FLAGS="-Wno-error=incompatible-pointer-types" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
    -DNOFORTRAN=ON \
&& cmake --build build --parallel $CORES && cmake --install build && cp -r install/* .. && cd ../..

# compile libtorch
cd external/pytorch-v2.5.0 && cmake -B build \
    -DBUILD_PYTHON=OFF \
    -DBUILD_SHARED_LIBS=$SHARED \
    -DCMAKE_BUILD_TYPE=${TYPE^} \
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
&& cmake --build build --parallel $CORES && cmake --install build && cp -r build/lib .. && cp -r install/* .. && cd ../..

# remove libraries of the other type
[ $SHARED == 1 ] && rm -f external/lib/*.a*
[ $STATIC == 1 ] && rm -f external/lib/*.s*

# remove sources
cd external && rm -rf eigen-3.4.0 fftw-3.3.10 libint-2.9.0 numactl-2.0.18 OpenBLAS-0.3.28 pytorch-v2.5.0 bin share *.tar.gz ; cd ..
