#!/bin/bash

CORES=$(nproc --all); TARGETS=("x86_64-linux")

export C_INCLUDE_PATH="$PWD/llvm/install/include:$C_INCLUDE_PATH"; export CMAKE_PREFIX_PATH="$PWD/eigen/install:$CMAKE_PREFIX_PATH"; export PATH="$PWD/zig-bin:$PATH"

for TARGET in "${TARGETS[@]}"; do

    echo "zig cc --target=$TARGET-gnu \"\$@\"" > zigcc && chmod +x zigcc && echo "zig c++ --target=$TARGET-gnu \"\$@\"" > zigcpp && chmod +x zigcpp

    tar -xzf    boost.tar.gz && mv boost-1.88.0               boost
    tar -xzf    eigen.tar.gz && mv eigen-3.4.0                eigen
    tar -xzf   libint.tar.gz && mv libint-2.11.1             libint
    tar -xzf     llvm.tar.gz && mv llvm-project-llvmorg-20.1.8 llvm
    tar -xzf openblas.tar.gz && mv OpenBLAS-0.3.30         openblas

    cd llvm && mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$PWD/../../zigcc" -DCMAKE_CXX_COMPILER="$PWD/../../zigcpp" -DCMAKE_INSTALL_PREFIX="$PWD/../install" -DLIBOMP_OMPD_SUPPORT=False -DLIBOMP_ENABLE_SHARED=False ../openmp && cmake --build . --parallel $CORES --verbose && cmake --install . && cd ../..

    cd openblas && make CC="$PWD/../zigcc" DYNAMIC_ARCH=1 HOSTCC=gcc NO_SHARED=1 NOFORTRAN=1 NUM_THREADS=128 USE_OPENMP=1 -j $CORES libs shared && make NO_SHARED=1 PREFIX="$PWD/install" install && cd ..

    cd boost && mkdir -p install && mkdir -p install/include && ./bootstrap.sh && ./b2 tools/bcp && ./dist/bin/bcp preprocessor install/include && rm -rf install/include/libs && cd ..

    cd eigen && mkdir -p build && cd build && cmake -DCMAKE_INSTALL_PREFIX="$PWD/../install" .. && cmake --build . --verbose && cmake --install . && cd ../install && mv include/eigen3/* include && rm -rf include/eigen3 && cd ../..

    cd libint && cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="$PWD/../zigcpp" -DCMAKE_INSTALL_PREFIX="$PWD/install" . && cmake --build . --parallel $CORES --verbose && cmake --install . && cd ..

    mkdir external-$TARGET && cp -r boost/install/* eigen/install/* libint/install/* llvm/install/* openblas/install/* external-$TARGET

    cd external-$TARGET && rm -rf bin share lib/*archer* lib/*cblas* lib/*la lib/cmake lib/pkgconfig && cd ..

    rm -rf boost eigen libint llvm openblas zigcc zigcpp

    wget -q -O external-$TARGET/include/exprtk.hpp https://raw.githubusercontent.com/ArashPartow/exprtk/refs/heads/master/exprtk.hpp

done
