#!/bin/bash

CORES=$(nproc --all)

export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$PWD/eigen/install"; export PATH="$PWD/bin:$PATH"

echo 'zig cc  "$@"' > zigcc  && chmod +x zigcc
echo 'zig c++ "$@"' > zigcpp && chmod +x zigcpp

wget -q --show-progress -O    eigen.tar.gz https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
wget -q --show-progress -O     fftw.tar.gz https://www.fftw.org/fftw-3.3.10.tar.gz
wget -q --show-progress -O      gsl.tar.gz https://sunsite.icm.edu.pl/pub/gnu/gsl/gsl-2.8.tar.gz
wget -q --show-progress -O   libint.tar.gz https://github.com/evaleev/libint/releases/download/v2.10.2/libint-2.10.2.tgz
wget -q --show-progress -O openblas.tar.gz https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.29/OpenBLAS-0.3.29.tar.gz

tar -xzf eigen.tar.gz    && rm    eigen.tar.gz
tar -xzf fftw.tar.gz     && rm     fftw.tar.gz
tar -xzf gsl.tar.gz      && rm      gsl.tar.gz
tar -xzf libint.tar.gz   && rm   libint.tar.gz
tar -xzf openblas.tar.gz && rm openblas.tar.gz

rm -rf eigen fftw gsl libint openblas && mv eigen* eigen && mv fftw* fftw ; mv gsl* gsl ; mv libint* libint ; mv OpenBLAS* openblas

cd fftw     && ./configure CC="$PWD/../zigcc" --enable-shared --prefix="$PWD/install" && make                               -j $CORES && make                       install && cd ..
cd gsl      && ./configure CC="$PWD/../zigcc"                 --prefix="$PWD/install" && make                               -j $CORES && make                       install && cd ..
cd openblas &&                                                                           make CC="$PWD/../zigcc" HOSTCC=gcc -j $CORES && make PREFIX="$PWD/install" install && cd ..

cd eigen && mkdir -p build && cd build && cmake -DCMAKE_INSTALL_PREFIX="$PWD/../install" .. && cmake --build . --parallel $CORES --verbose && cmake --install . && cd ../..

cd libint && cmake -DBUILD_SHARED_LIBS=ON  -DCMAKE_CXX_COMPILER="$PWD/../zigcpp" -DCMAKE_INSTALL_PREFIX="$PWD/install" . && cmake --build . --parallel $CORES --verbose && cmake --install . && cd ..
cd libint && cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_CXX_COMPILER="$PWD/../zigcpp" -DCMAKE_INSTALL_PREFIX="$PWD/install" . && cmake --build . --parallel $CORES --verbose && cmake --install . && cd ..

rm -rf external && mkdir external && cp -r eigen/install/* fftw/install/* gsl/install/* libint/install/* openblas/install/* external

mv external/include/eigen3/Eigen external/include && rm -rf eigen fftw gsl libint openblas zigcc zigcpp external/include/eigen3
