#!/bin/bash

# download library
mkdir -p external && git clone https://github.com/boostorg/boost external/libboost

# checkout library
cd external/libboost && git checkout a7090e8ce184501cfc9e80afa6cafb5bfd3b371c && cd -

# download submodules
cd external/libboost && git submodule update --init && cd -

# configure library
cd external/libboost && ./bootstrap.sh --prefix="$PWD/install" --with-libraries="atomic" && cd -

# compile and install library
cd external/libboost && ./b2 install && cd -

# copy the compiled library
cp -r external/libboost/install/include external/libboost/install/lib external/

# remove redundant files
# rm -rf external/libboost external/lib/cmake external/lib/*.so*
