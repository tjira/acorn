#!/bin/bash

# clone the repository
git clone --depth 1 https://github.com/evaleev/libint.git libint

# checkout to the specific commit
cd libint && git checkout fe4ab095dfbc1fb84b18e51ee11f0f0071b5905e && cd -

# compile libint
cd libint && ./autogen.sh && ./configure CXX=g++ CXXFLAGS="-mavx -s -O3" --prefix="$PWD/install" --enable-1body=1 --enable-eri=1 --with-max-am=6 && make -j2 && make install && cd ..
