#!/bin/bash

CORES=$(nproc --all); TARGETS=("x86_64-linux")

export C_INCLUDE_PATH="$C_INCLUDE_PATH:$PWD/llvm/install/include"; export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$PWD/eigen/install"; export PATH="$PWD/zig-bin:$PATH"

for TARGET in "${TARGETS[@]}"; do

    rm -rf zig-bin && mkdir zig-bin && ZIG_VERSION="0.14.1" && ZLS_VERSION="0.14.0"

    wget -q -O zig.tar.xz https://ziglang.org/download/$ZIG_VERSION/zig-x86_64-linux-$ZIG_VERSION.tar.xz
    wget -q -O zls.tar.xz https://github.com/zigtools/zls/releases/download/$ZLS_VERSION/zls-x86_64-linux.tar.xz

    tar -xf zig.tar.xz && mv zig-x86_64-linux*/lib zig-x86_64-linux*/zig zig-bin && tar -xf zls.tar.xz -C zig-bin && rm -rf zig-x86_64-linux* *.tar.xz

    echo "zig cc --target=$TARGET-gnu \"\$@\"" > zigcc && chmod +x zigcc && echo "zig c++ --target=$TARGET-gnu \"\$@\"" > zigcpp && chmod +x zigcpp

    wget -q -O boost.tar.gz    https://github.com/boostorg/boost/releases/download/boost-1.88.0/boost-1.88.0-b2-nodocs.tar.gz
    wget -q -O eigen.tar.gz    https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    wget -q -O libint.tar.gz   https://github.com/evaleev/libint/releases/download/v2.10.2/libint-2.10.2.tgz
    wget -q -O llvm.tar.gz     https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-20.1.4.tar.gz
    wget -q -O openblas.tar.gz https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.29/OpenBLAS-0.3.29.tar.gz

    tar -xzf boost.tar.gz    && rm    boost.tar.gz
    tar -xzf eigen.tar.gz    && rm    eigen.tar.gz
    tar -xzf libint.tar.gz   && rm   libint.tar.gz
    tar -xzf llvm.tar.gz     && rm     llvm.tar.gz
    tar -xzf openblas.tar.gz && rm openblas.tar.gz

    rm -rf boost eigen libint llvm openblas && mv boost* boost ; mv eigen* eigen ; mv libint* libint ; mv llvm* llvm ; mv OpenBLAS* openblas

    cd llvm && mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$PWD/../../zigcc" -DCMAKE_CXX_COMPILER="$PWD/../../zigcpp" -DCMAKE_INSTALL_PREFIX="$PWD/../install" -DLIBOMP_OMPD_SUPPORT=False -DLIBOMP_ENABLE_SHARED=False ../openmp && cmake --build . --parallel $CORES --verbose && cmake --install . && cd ../..

    cd openblas && make CC="$PWD/../zigcc" DYNAMIC_ARCH=1 HOSTCC=gcc NO_SHARED=1 NOFORTRAN=1 NUM_THREADS=128 USE_OPENMP=1 -j $CORES libs shared && make NO_SHARED=1 PREFIX="$PWD/install" install && cd ..

    cd boost && mkdir -p install && mkdir -p install/include && ./bootstrap.sh && ./b2 tools/bcp && ./dist/bin/bcp preprocessor install/include && rm -rf install/include/libs && cd ..

    cd eigen && mkdir -p build && cd build && cmake -DCMAKE_INSTALL_PREFIX="$PWD/../install" .. && cmake --build . --verbose && cmake --install . && cd ../install && mv include/eigen3/* include && rm -rf include/eigen3 && cd ../..

    cd libint && cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="$PWD/../zigcpp" -DCMAKE_INSTALL_PREFIX="$PWD/install" . && cmake --build . --parallel $CORES --verbose && cmake --install . && cd ..

    rm -rf external-$TARGET && mkdir external-$TARGET && cp -r boost/install/* eigen/install/* libint/install/* llvm/install/* openblas/install/* external-$TARGET

    cd external-$TARGET && rm -rf bin share lib/*archer* lib/*cblas* lib/*la lib/cmake lib/pkgconfig && cd ..

    rm -rf boost eigen libint llvm openblas zigcc zigcpp

    wget -q -O external-$TARGET/include/exprtk.hpp https://raw.githubusercontent.com/ArashPartow/exprtk/refs/heads/master/exprtk.hpp

done
