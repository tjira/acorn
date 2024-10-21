#!/bin/bash

# compile parameters
SHARED="OFF"; CORES=2; if [ $# -ne 0 ] && [ $# -ne 2 ]; then echo "ARGUMENTS: SHARED(ON/OFF) CORES(1,2,...)"; exit 1; fi; if [ $# -eq 2 ]; then SHARED=$1; CORES=$2; fi

# download the library
mkdir -p external && wget -O external/libfftw.tar.gz https://www.fftw.org/fftw-3.3.10.tar.gz && cd external && rm -rf libfftw && tar -xzvf libfftw* && mv fftw* libfftw && cd -

# configure the library
cd external/libfftw && cmake -B build -DBUILD_SHARED_LIBS=$SHARED -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install" && cd -

# compile and install the library
cd external/libfftw && cmake --build build --parallel $CORES && cmake --install build && cd -

# copy the compiled library
cp -r external/libfftw/install/include external/libfftw/install/lib external/
