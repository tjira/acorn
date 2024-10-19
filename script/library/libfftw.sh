#!/bin/bash

# download the library
mkdir -p external && wget -O external/libfftw.tar.gz https://www.fftw.org/fftw-3.3.10.tar.gz && cd external && rm -rf libfftw && tar -xzvf libfftw* && mv fftw* libfftw && cd -

# configure the library
cd external/libfftw && cmake -B build -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX="$PWD/install" && cd -

# compile and install the library
cd external/libfftw && cmake --build build --parallel 2 && cmake --install build && cd -

# copy the compiled library
cp -r external/libfftw/install/include external/libfftw/install/lib external/

# rename the compiled library
cd external/lib && mv libfftw3.a libfftw.a ; cd -
