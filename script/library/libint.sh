#!/bin/bash

# compile parameters
SHARED="no"; STATIC="yes" CORES=2; if [ $# -ne 0 ] && [ $# -ne 2 ]; then echo "ARGUMENTS: SHARED(ON/OFF) CORES(1,2,...)"; exit 1; fi; if [ $# -eq 2 ]; then SHARED=$1; CORES=$2; fi; if [ $SHARED == "ON" ]; then SHARED="yes"; STATIC="no"; fi

# download the library
mkdir -p external && git clone https://github.com/evaleev/libint.git external/libint

# export the paths to dependencies
export CPLUS_INCLUDE_PATH="$PWD/external/include:$CPLUS_INCLUDE_PATH"; export LIBRARY_PATH="$PWD/external/lib:$LIBRARY_PATH"

# configure the library
cd external/libint && ./autogen.sh && ./configure CXX=c++ CXXFLAGS="-O3 -march=native" --prefix="$PWD/install" --enable-1body=1 --enable-eri=1 --enable-shared=$SHARED --enable-static=$STATIC --with-cxxgen-optflags="-O3 -march=native" --with-max-am=6 && cd -

# compile and install the library
cd external/libint && make -j$CORES && make install && cd -

# copy the compiled library
cp -r external/libint/install/include external/libint/install/lib external/
