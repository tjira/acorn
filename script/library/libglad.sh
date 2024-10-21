#!/bin/bash

# download the library
mkdir -p external && git clone https://github.com/Dav1dde/glad.git external/libglad

# configure the library
cd external/libglad && pip install -r requirements.txt && python -m glad --api="gl:core=4.2" --out-path="$PWD/install" && cd -

# compile and install the library
cd external/libglad && mkdir -p install/lib && gcc -c install/src/gl.c -o install/lib/libglad.a -Iinstall/include && cd -

# copy the compiled library
cp -r external/libglad/install/include external/libglad/install/lib external/

# remove the source
rm -rf external/libglad*
