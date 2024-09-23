#!/bin/bash

# download libint
mkdir -p external && git clone --recursive https://github.com/boostorg/boost external/libboost

# configure boost
cd external/libboost && ./bootstrap.sh --prefix="$PWD/install" --with-libraries="atomic" && cd -

# compile and install boost
cd external/libboost && ./b2 install && cd -

# copy the compiled library
cp -r external/libboost/install/include external/libboost/install/lib external/

# remove redundant files
rm -rf external/libboost external/lib/cmake external/lib/*.so*
