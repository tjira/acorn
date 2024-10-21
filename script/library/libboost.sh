#!/bin/bash

# download the library
mkdir -p external && git clone https://github.com/boostorg/boost external/libboost && cd external/libboost && git submodule update --init && cd -

# configure the library
cd external/libboost && ./bootstrap.sh --prefix="$PWD/install" --with-libraries="atomic" && cd -

# compile and install the library
cd external/libboost && ./b2 install && cd -

# copy the compiled library
cp -r external/libboost/install/include external/libboost/install/lib external/

# remove the source
rm -rf external/libboost*
