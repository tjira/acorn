#!/bin/bash

# compile parameters
SHARED="OFF"; CORES=2; if [ $# -ne 0 ] && [ $# -ne 2 ]; then echo "ARGUMENTS: SHARED(ON/OFF) CORES(1,2,...)"; exit 1; fi; if [ $# -eq 2 ]; then SHARED=$1; CORES=$2; fi

# download the library
mkdir -p external && git clone https://github.com/g-truc/glm.git external/libglm

# configure the library
cd external/libglm && cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=$SHARED -DCMAKE_INSTALL_PREFIX="$PWD/install" && cd -

# compile and install the library
cd external/libglm && cmake --build build --parallel $CORES && cmake --install build && cd -

# copy the compiled library
cp -r external/libglm/install/include external/libglm/install/lib external/
