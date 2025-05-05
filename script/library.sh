#!/bin/bash

CORES=$(nproc --all)

export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$PWD/eigen/install"; export PATH="$PWD/bin:$PATH"

echo 'zig cc  "$@"' > zigcc  && chmod +x zigcc && echo 'zig c++ "$@"' > zigcpp && chmod +x zigcpp

wget -q -O    boost.tar.gz https://archives.boost.io/release/1.88.0/source/boost_1_88_0.tar.gz
wget -q -O    eigen.tar.gz https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
wget -q -O     fftw.tar.gz https://www.fftw.org/fftw-3.3.10.tar.gz
wget -q -O      gsl.tar.gz https://sunsite.icm.edu.pl/pub/gnu/gsl/gsl-2.8.tar.gz
wget -q -O   libint.tar.gz https://github.com/evaleev/libint/releases/download/v2.10.2/libint-2.10.2.tgz
wget -q -O openblas.tar.gz https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.29/OpenBLAS-0.3.29.tar.gz
wget -q -O     llvm.tar.gz https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-20.1.4.tar.gz

tar -xzf llvm.tar.gz     && rm     llvm.tar.gz
tar -xzf boost.tar.gz    && rm    boost.tar.gz
tar -xzf eigen.tar.gz    && rm    eigen.tar.gz
tar -xzf fftw.tar.gz     && rm     fftw.tar.gz
tar -xzf gsl.tar.gz      && rm      gsl.tar.gz
tar -xzf libint.tar.gz   && rm   libint.tar.gz
tar -xzf openblas.tar.gz && rm openblas.tar.gz

rm -rf boost eigen fftw gsl libint llvm openblas && mv boost* boost ; mv eigen* eigen ; mv fftw* fftw ; mv gsl* gsl ; mv libint* libint ; mv llvm* llvm ; mv OpenBLAS* openblas

cd llvm && mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/../install" -DLIBOMP_ENABLE_SHARED=True  ../openmp && cmake --build . --parallel $CORES --verbose && cmake --install . && cd ../..
cd llvm && mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/../install" -DLIBOMP_ENABLE_SHARED=False ../openmp && cmake --build . --parallel $CORES --verbose && cmake --install . && cd ../..

cd boost    && ./bootstrap.sh --prefix="$PWD/install" --with-libraries="atomic" && ./b2 install && cd ..

cd fftw     && ./configure --enable-openmp --enable-shared --prefix="$PWD/install" && cd ..
cd gsl      && ./configure                                 --prefix="$PWD/install" && cd ..

cd fftw     && make                              -j $CORES && make                       install && cd ..
cd gsl      && make                              -j $CORES && make                       install && cd ..
cd openblas && make NUM_THREADS=128 USE_OPENMP=1 -j $CORES && make PREFIX="$PWD/install" install && cd ..

cd eigen && mkdir -p build && cd build && cmake -DCMAKE_INSTALL_PREFIX="$PWD/../install" .. && cmake --build . --parallel $CORES --verbose && cmake --install . && cd ../..

cd libint && cmake -DBUILD_SHARED_LIBS=ON  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="$PWD/../zigcpp" -DCMAKE_INSTALL_PREFIX="$PWD/install" . && cmake --build . --parallel $CORES --verbose && cmake --install . && cd ..
cd libint && cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="$PWD/../zigcpp" -DCMAKE_INSTALL_PREFIX="$PWD/install" . && cmake --build . --parallel $CORES --verbose && cmake --install . && cd ..

rm -rf external && mkdir external && cp -r boost/install/* eigen/install/* fftw/install/* gsl/install/* libint/install/* llvm/install/* openblas/install/* external

mv external/include/eigen3/Eigen external/include/eigen3/unsupported external/include && rm -rf boost eigen fftw gsl libint llvm openblas zigcc zigcpp external/include/eigen3
