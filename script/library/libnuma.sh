#!/bin/bash

# compile parameters
SHARED="no"; STATIC="yes" CORES=2; if [ $# -ne 0 ] && [ $# -ne 2 ]; then echo "ARGUMENTS: SHARED(ON/OFF) CORES(1,2,...)"; exit 1; fi; if [ $# -eq 2 ]; then SHARED=$1; CORES=$2; fi; if [ $SHARED == "ON" ]; then SHARED="yes"; STATIC="no"; fi

# download the library
mkdir -p external && git clone https://github.com/numactl/numactl.git external/libnuma

# configure the library
cd external/libnuma && ./autogen.sh && ./configure --prefix="$PWD/install" --enable-shared=$SHARED --enable-static=$STATIC && cd -

# compile and install the library
cd external/libnuma && make -j$CORES && make install && cd -

# copy the compiled library
cp -r external/libnuma/install/include external/libnuma/install/lib external/
