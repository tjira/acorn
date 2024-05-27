#!/bin/bash

# download and extract FFTW
wget https://www.fftw.org/fftw-3.3.10.tar.gz && tar -xzvf fftw-3.3.10.tar.gz && rm fftw-3.3.10.tar.gz && mv fftw-3.3.10 libfftw

# configure FFTW
cd libfftw && cmake -B build -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTS=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX="$PWD/install" && cd ..

# compile and install FFTW
cd libfftw && cmake --build build && cmake --install build && cd ..
