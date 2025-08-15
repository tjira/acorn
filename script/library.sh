#!/bin/bash

CORES=$(nproc --all); DIRNAME=$(dirname "$(realpath $0)"); TARGETS=("x86_64-linux")

CMAKE_GENERAL=(
    -DBUILD_SHARED_LIBS=OFF
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_C_COMPILER="$DIRNAME/../zigcc"
    -DCMAKE_CXX_COMPILER="$DIRNAME/../zigcpp"
)

CMAKE_LIBINT=(
    -DLIBINT_USE_BUNDLED_BOOST=True
)

CMAKE_OMP=(
    -DLIBOMP_ENABLE_SHARED=False
    -DLIBOMP_OMPD_SUPPORT=False
    -DLIBOMP_OMPT_SUPPORT=False
)

MAKE_OPENBLAS=(
    AR="$DIRNAME/../zigar"
    CC="$DIRNAME/../zigcc"
    DYNAMIC_ARCH=1
    HOSTCC=gcc
    NO_SHARED=1
    NOFORTRAN=1
    NUM_THREADS=128
    RANLIB="$DIRNAME/../zigranlib"
    USE_OPENMP=1
)

cmake_install() {
    cmake --build . --parallel $CORES --verbose && cmake --install .
}

for TARGET in "${TARGETS[@]}"; do

    export C_INCLUDE_PATH="$DIRNAME/../external-$TARGET/include:$C_INCLUDE_PATH"

    CMAKE_GENERAL+=(
        -DCMAKE_INSTALL_PREFIX="$DIRNAME/../external-$TARGET"
        -DCMAKE_PREFIX_PATH="$DIRNAME/../external-$TARGET"
    )

    echo "$DIRNAME/../zig-bin/zig ar                          \"\$@\"" > zigar     && chmod +x     zigar
    echo "$DIRNAME/../zig-bin/zig cc     --target=$TARGET-gnu \"\$@\"" > zigcc     && chmod +x     zigcc
    echo "$DIRNAME/../zig-bin/zig c++    --target=$TARGET-gnu \"\$@\"" > zigcpp    && chmod +x    zigcpp
    echo "$DIRNAME/../zig-bin/zig ranlib                      \"\$@\"" > zigranlib && chmod +x zigranlib

    tar -xzf    external-$TARGET/eigen.tar.gz && mv eigen*       eigen
    tar -xzf   external-$TARGET/libint.tar.gz && mv libint*     libint
    tar -xf      external-$TARGET/llvm.tar.xz && mv llvm*         llvm
    tar -xzf external-$TARGET/openblas.tar.gz && mv OpenBLAS* openblas

    cd llvm   && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}" "${CMAKE_OMP[@]}"    ../openmp && cmake_install && cd ../..
    cd eigen  && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}"                      ..        && cmake_install && cd ../..
    cd libint && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}" "${CMAKE_LIBINT[@]}" ..        && cmake_install && cd ../..

    cd openblas && make "${MAKE_OPENBLAS[@]}" -j $CORES libs shared && make NO_SHARED=1 PREFIX="$DIRNAME/../external-$TARGET" install && cd ..

    rm -rf eigen libint llvm openblas zigar zigcc zigcpp zigranlib

done
