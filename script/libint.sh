#!/bin/bash

# clone the repository
git clone --depth 1 https://github.com/evaleev/libint.git libint

# compile libint
cd libint && ./autogen.sh && ./configure --enable-1body=1 --enable-eri=1 --prefix="$PWD/install" --with-max-am=4 && make -j2 && make install && cd ..
