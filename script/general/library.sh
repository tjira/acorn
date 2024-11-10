#!/bin/bash

# check the number of arguments
[ $# -ne 2 ] && echo "USAGE: $0 MODE CORES" && exit 1

# assign the arguments
CORES=$2; MODE=$1; [ $MODE == "STATIC" ] && STATIC=1 || STATIC=0; [ $MODE == "SHARED" ] && SHARED=1 || SHARED=0; [ $SHARED -eq 0 ] && [ $STATIC -eq 0 ] && echo "INVALID MODE" && exit 1

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

# download glfw
wget -q -O external/glfw.zip https://github.com/glfw/glfw/releases/download/3.4/glfw-3.4.zip

# download glm
wget -q -O external/glm.zip https://github.com/g-truc/glm/releases/download/1.0.1/glm-1.0.1-light.zip

# download imgui
wget -q -O external/imgui.tar.gz https://github.com/ocornut/imgui/archive/refs/tags/v1.91.4.tar.gz

# download imguifiledialog
wget -q -O external/imguifiledialog.tar.gz https://github.com/aiekick/ImGuiFileDialog/archive/refs/tags/v0.6.7.tar.gz

# download implot
wget -q -O external/implot.tar.gz https://github.com/epezent/implot/archive/refs/tags/v0.16.tar.gz

# download libint
wget -q -O external/libint.tar.gz https://github.com/evaleev/libint/releases/download/v2.9.0/libint-2.9.0-mpqc4.tgz

# download libnuma
wget -q -O external/libnuma.tar.gz https://github.com/numactl/numactl/releases/download/v2.0.18/numactl-2.0.18.tar.gz

# download libopenblas
wget -q -O external/openblas.tar.gz https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.28/OpenBLAS-0.3.28.tar.gz

# download libtorch
wget -q -O external/libtorch.tar.gz https://github.com/pytorch/pytorch/releases/download/v2.5.0/pytorch-v2.5.0.tar.gz

# unpack the archives
cd external && for ARCHIVE in *.tar.gz; do tar -xzf  $ARCHIVE --warning=no-unknown-keyword; done; cd ..
cd external && for ARCHIVE in *.zip;    do unzip -qq $ARCHIVE                             ; done; cd ..

# compile eigen
cd external/eigen-3.4.0 && cmake -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
    -DEIGEN_BUILD_DOC=OFF \
    -DEIGEN_TEST_NOQT=ON \
&& cmake --build build --parallel $CORES && cmake --install build && cp -r install/* .. && cd ../..

# compile fftw
cd external/fftw-3.3.10 && cmake -B build \
    -DBUILD_SHARED_LIBS=$SHARED \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
&& cmake --build build --parallel $CORES && cmake --install build && cp -r install/* .. && cd ../..

# compile glfw
cd external/glfw-3.4 && cmake -B build \
    -DBUILD_SHARED_LIBS=$SHARED \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
&& cmake --build build --parallel $CORES && cmake --install build && cp -r install/* .. && cd ../..

# compile libint
cd external/libint-2.9.0 && cmake -B build \
    -DBUILD_SHARED_LIBS=$SHARED \
    -DCMAKE_BUILD_TYPE=Release \
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
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_FLAGS="-Wno-error=incompatible-pointer-types" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
    -DNOFORTRAN=ON \
&& cmake --build build --parallel $CORES && cmake --install build && cp -r install/* .. && cd ../..

# compile static openblas
[ $STATIC == 1 ] && cd external/OpenBLAS-0.3.28 && cmake -B build \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_FLAGS="-Wno-error=incompatible-pointer-types" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
    -DNOFORTRAN=ON \
&& cmake --build build --parallel $CORES && cmake --install build && cp -r install/* .. && cd ../..

# compile libtorch
cd external/pytorch-v2.5.0 && cmake -B build \
    -DBUILD_PYTHON=OFF \
    -DBUILD_SHARED_LIBS=$SHARED \
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
&& cmake --build build --parallel $CORES && cmake --install build && cp -r build/lib .. && cp -r install/* .. && cd ../..

# download glad
pip install glad && glad --api="gl=4.2" --generator="c" --profile="core" --out-path="$PWD/external"

# copy headers
cp -r external/imgui-1.91.4/*.h external/imgui-1.91.4/backends external/ImGuiFileDialog-0.6.7/*.h external/ImGuiFileDialog-0.6.7/dirent external/implot-0.16/*.h external/include && mv external/glm external/include

# copy sources
cp -r external/imgui-1.91.4/*.cpp external/imgui-1.91.4/backends external/ImGuiFileDialog-0.6.7/*.cpp external/implot-0.16/*.cpp external/src 

# remove libraries of the other type
[ $SHARED == 1 ] && rm -f external/lib/*.a*
[ $STATIC == 1 ] && rm -f external/lib/*.s*

# remove sources
cd external && rm -rf eigen-3.4.0 fftw-3.3.10 glfw-3.4 imgui-1.91.4 ImGuiFileDialog-0.6.7 implot-0.16 libint-2.9.0 numactl-2.0.18 OpenBLAS-0.3.28 pytorch-v2.5.0 bin share *.tar.gz *.zip ; cd ..
