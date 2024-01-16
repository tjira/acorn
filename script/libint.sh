#!/bin/bash

# clone the repository
git clone --branch v2.8.1 https://github.com/evaleev/libint.git libint

# compile libint
cd libint && ./autogen.sh && ./configure --prefix="$PWD/install" --enable-1body=1 --enable-eri=1 --with-max-am=6 --with-cxxgen-optflags="-mavx -s -O3" && make && make install && cd ..
