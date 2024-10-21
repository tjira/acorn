#!/bin/bash

# compile parameters
SHARED="no"; STATIC="yes" CORES=2; if [ $# -ne 0 ] && [ $# -ne 2 ]; then echo "ARGUMENTS: SHARED(ON/OFF) CORES(1,2,...)"; exit 1; fi; if [ $# -eq 2 ]; then SHARED=$1; CORES=$2; fi; if [ $SHARED == "ON" ]; then SHARED="yes"; STATIC="no"; fi

# download the library
mkdir -p external && wget -O external/libgmp.tar.xz https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz && cd external && rm -rf libgmp && tar -xvf libgmp* && mv gmp* libgmp && cd -

# configure the library
cd external/libgmp && ./configure --enable-cxx --prefix="$PWD/install" --enable-shared=$SHARED --enable-static=$STATIC && cd -

# compile and install the library
cd external/libgmp && make -j$CORES && make install && cd -

# copy the compiled library
cp -r external/libgmp/install/include external/libgmp/install/lib external/

# remove the source
rm -rf external/libgmp*
