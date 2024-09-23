#!/bin/bash

# download GMP
mkdir -p external && wget -O external/libgmp.tar.xz https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz

# unpack GMP
cd external && rm -rf libgmp && tar -xvf libgmp.tar.xz && mv gmp* libgmp && cd -

# configure GMP
cd external/libgmp && ./configure --enable-cxx --prefix="$PWD/install" && cd -

# compile and install GMP
cd external/libgmp && make -j2 && make install && cd -

# copy the compiled library
cp -r external/libgmp/install/include external/libgmp/install/lib external/

# remove redundant files
rm -rf external/libgmp external/lib/pkgconfig external/lib/*.la external/lib/*.so*
