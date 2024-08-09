#!/bin/bash

# download FFTW
mkdir -p external && wget -O external/libfftw.tar.gz https://www.fftw.org/fftw-3.3.10.tar.gz

# unpack FFTW
cd external && rm -rf libfftw && tar -xzvf libfftw.tar.gz && mv fftw* libfftw && cd -

# configure FFTW
cd external/libfftw && cmake -B build -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install" && cd -

# compile and install FFTW
cd external/libfftw && cmake --build build --parallel 2 && cmake --install build && cd -

# copy the compiled library
cp -r external/libfftw/install/include external/libfftw/install/lib external/

# remove redundant files
rm -rf external/libfftw external/include/*.f external/include/*.f03 external/lib/cmake external/lib/pkgconfig

# rename the shared library
cd external/lib && [[ -f libfftw3.so.3 ]] && mv libfftw3.so.3 libfftw.so && patchelf --set-soname libfftw.so libfftw.so && rm libfftw3.* ; cd -

# rename the static library
cd external/lib && [[ -f libfftw3.a ]] && mv libfftw3.a libfftw.a ; cd -
