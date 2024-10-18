#!/bin/bash

# download the library
mkdir -p external && wget -O external/libgmp.tar.xz https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz && cd external && rm -rf libgmp && tar -xvf libgmp* && mv gmp* libgmp && cd -

# configure the library
cd external/libgmp && ./configure --enable-cxx --prefix="$PWD/install" && cd -

# compile and install the library
cd external/libgmp && make -j2 && make install && cd -

# copy the compiled library
cp -r external/libgmp/install/include external/libgmp/install/lib external/
