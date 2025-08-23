#!/bin/bash

if [ $# -lt 1 ]; then exit 1; fi

CORES=$(nproc --all); DIRNAME=$(dirname "$(realpath $0)"); TARGET=${1,,}; HOST="$(uname -m)-$(uname -s | tr '[:upper:]' '[:lower:]')"

CMAKE_GENERAL=(
    -DBUILD_SHARED_LIBS=False
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX="$PWD/external-$TARGET"
    -DCMAKE_PREFIX_PATH="$PWD/external-$TARGET"
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
    HOSTCC=gcc
    NOFORTRAN=1
    NO_SHARED=1
    NUM_THREADS=128
    PREFIX="$PWD/external-$TARGET"
    USE_OPENMP=1
)

cmake_install() {
    cmake --build . --parallel $CORES --verbose && cmake --install .
}

if [[ "$HOST" == "$TARGET" ]]; then MAKE_OPENBLAS+=(DYNAMIC_ARCH=1); else MAKE_OPENBLAS+=(CROSS=1 TARGET=GENERIC); fi

export C_INCLUDE_PATH="$PWD/external-$TARGET/include:$C_INCLUDE_PATH"

echo -e "#!/usr/bin/env bash\n\n$PWD/zig-bin/zig ar                          \"\$@\"" > zigar     && chmod +x     zigar
echo -e "#!/usr/bin/env bash\n\n$PWD/zig-bin/zig cc     --target=$TARGET-gnu \"\$@\"" > zigcc     && chmod +x     zigcc
echo -e "#!/usr/bin/env bash\n\n$PWD/zig-bin/zig c++    --target=$TARGET-gnu \"\$@\"" > zigcpp    && chmod +x    zigcpp
echo -e "#!/usr/bin/env bash\n\n$PWD/zig-bin/zig ranlib                      \"\$@\"" > zigranlib && chmod +x zigranlib

export CC="$PWD/zigcc"; export CXX="$PWD/zigcpp"; export AR="$PWD/zigar"; export RANLIB="$PWD/zigranlib"

rm -rf eigen    && tar -xzf    external-$HOST/eigen.tar.gz && mv eigen*       eigen
rm -rf libint   && tar -xzf   external-$HOST/libint.tar.gz && mv libint*     libint
rm -rf llvm     && tar -xf      external-$HOST/llvm.tar.xz && mv llvm*         llvm
rm -rf openblas && tar -xzf external-$HOST/openblas.tar.gz && mv OpenBLAS* openblas

cd llvm   && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}" "${CMAKE_OMP[@]}"    ../openmp && cmake_install && cd ../..
cd eigen  && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}"                      ..        && cmake_install && cd ../..
cd libint && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}" "${CMAKE_LIBINT[@]}" ..        && cmake_install && cd ../..

cd openblas && make "${MAKE_OPENBLAS[@]}" -j $CORES libs shared && make "${MAKE_OPENBLAS[@]}" install && cd ..

wget -q -O external-$TARGET/include/exprtk.hpp https://raw.githubusercontent.com/ArashPartow/exprtk/cc1b800c2bd1ac3ac260478c915d2aec6f4eb41c/exprtk.hpp

rm -rf eigen libint llvm openblas zigar zigcc zigcpp zigranlib
