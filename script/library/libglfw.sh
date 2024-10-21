#!/bin/bash

# compile parameters
SHARED="OFF"; CORES=2; if [ $# -ne 0 ] && [ $# -ne 2 ]; then echo "ARGUMENTS: SHARED(ON/OFF) CORES(1,2,...)"; exit 1; fi; if [ $# -eq 2 ]; then SHARED=$1; CORES=$2; fi

# download the library
mkdir -p external && git clone https://github.com/glfw/glfw.git external/libglfw

# configure the library
cd external/libglfw && cmake -B build -DBUILD_SHARED_LIBS=$SHARED -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install" && cd -

# compile and install the library
cd external/libglfw && cmake --build build --parallel $CORES && cmake --install build && cd -

# copy the compiled library
cp -r external/libglfw/install/include external/libglfw/install/lib external/

# remove the source
rm -rf external/libglfw*
