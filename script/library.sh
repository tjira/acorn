#!/bin/bash

CORES=$(nproc --all); DIRNAME=$(dirname "$(realpath $0)"); TARGETS=("x86_64-linux")

export C_INCLUDE_PATH="$PWD/llvm/install/include:$C_INCLUDE_PATH"; export CMAKE_PREFIX_PATH="$PWD/eigen/install:$CMAKE_PREFIX_PATH"; export PATH="$PWD/zig-bin:$PATH"

CMAKE_GENERAL=(
    -DBUILD_SHARED_LIBS=OFF
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_C_COMPILER="$DIRNAME/../zigcc"
    -DCMAKE_CXX_COMPILER="$DIRNAME/../zigcpp"
    -DCMAKE_INSTALL_PREFIX="../install"
)

CMAKE_LIBINT=(
    -DLIBINT_USE_BUNDLED_BOOST=True
)

CMAKE_OMP=(
    -DLIBOMP_ENABLE_SHARED=False
    -DLIBOMP_OMPD_SUPPORT=False
)

cmake_install() {
    cmake --build . --parallel $CORES --verbose && cmake --install .
}

for TARGET in "${TARGETS[@]}"; do

    echo "zig cc --target=$TARGET-gnu \"\$@\"" > zigcc && chmod +x zigcc && echo "zig c++ --target=$TARGET-gnu \"\$@\"" > zigcpp && chmod +x zigcpp

    tar -xzf    external-$TARGET/eigen.tar.gz && mv eigen*       eigen
    tar -xzf   external-$TARGET/libint.tar.gz && mv libint*     libint
    tar -xf      external-$TARGET/llvm.tar.xz && mv llvm*         llvm
    tar -xzf external-$TARGET/openblas.tar.gz && mv OpenBLAS* openblas

    cd llvm   && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}" "${CMAKE_OMP[@]}"    ../openmp && cmake_install && cd ../..
    cd eigen  && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}"                      ..        && cmake_install && cd ../..
    cd libint && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}" "${CMAKE_LIBINT[@]}" ..        && cmake_install && cd ../..

    cd openblas && make CC="$PWD/../zigcc" DYNAMIC_ARCH=1 HOSTCC=gcc NO_SHARED=1 NOFORTRAN=1 NUM_THREADS=128 USE_OPENMP=1 -j $CORES libs shared && make NO_SHARED=1 PREFIX="$PWD/install" install && cd ..

    cp -r libint/install/* llvm/install/* openblas/install/* external-$TARGET && cp -r eigen/install/include/eigen3/Eigen eigen/install/include/eigen3/unsupported external-$TARGET/include

    rm -rf eigen libint llvm openblas zigcc zigcpp && cd external-$TARGET && rm -rf bin share lib/*archer* lib/*cblas* lib/*la lib/cmake lib/pkgconfig && cd ..

done
