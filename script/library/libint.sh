#!/bin/bash

# download the library
mkdir -p external && git clone https://github.com/evaleev/libint.git external/libint

# export the paths to dependencies
export CPLUS_INCLUDE_PATH="$PWD/external/include:$CPLUS_INCLUDE_PATH"; export LIBRARY_PATH="$PWD/external/lib:$LIBRARY_PATH"

# configure the library
cd external/libint && ./autogen.sh && ./configure CXX=g++ CXXFLAGS="-march=native -s -O3 -flto=auto" --prefix="$PWD/install" --enable-1body=1 --enable-eri=1 --with-max-am=6 && cd -

# compile and install the library
cd external/libint && make -j2 && make install && cd -

# copy the compiled library
cp -r external/libint/install/include external/libint/install/lib external/

# rename the compiled library
cd external/lib && mv libint2.a libint.a ; cd -
