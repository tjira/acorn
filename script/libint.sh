#!/bin/bash

# download libint
mkdir -p external && wget -O external/libint.tgz https://github.com/evaleev/libint/releases/download/v2.9.0/libint-2.9.0-mpqc4.tgz

# unpack libint
cd external && rm -rf libint && tar -xvf libint.tgz --warning=no-unknown-keyword && mv libint-* libint && cd -

# configure libint
cd external/libint && cmake -B build -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install" && cd -

# compile and install libint
cd external/libint && cmake --build build --parallel 2 && cmake --install build && cd -

# copy the compiled library
cp -r external/libint/install/include external/libint/install/lib external/

# remove redundant files
rm -rf external/libint external/lib/cmake external/lib/pkgconfig external/lib/libint2.so.2

# move links
mv external/lib/libint2.so.2.9.0 external/lib/libint2.so
