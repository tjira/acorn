#!/bin/bash

# clone the repository
git clone https://github.com/evaleev/libint.git libint

# compile libint
cd libint && ./autogen.sh && ./configure CXX=g++ CXXFLAGS="-march=native -s -O3" --prefix="$PWD/install" --enable-1body=1 --enable-eri=1 --with-max-am=6 && make -j2 && make install && cd ..
