#!/bin/bash

wget -O     fftw.tar.gz https://www.fftw.org/fftw-3.3.10.tar.gz
wget -O      gsl.tar.gz https://sunsite.icm.edu.pl/pub/gnu/gsl/gsl-2.8.tar.gz
wget -O   libint.tar.gz https://github.com/evaleev/libint/releases/download/v2.10.2/libint-2.10.2.tgz
wget -O openblas.tar.gz https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.29/OpenBLAS-0.3.29.tar.gz

tar -xzvf fftw.tar.gz     && rm     fftw.tar.gz
tar -xzvf gsl.tar.gz      && rm      gsl.tar.gz
tar -xzvf libint.tar.gz   && rm   libint.tar.gz
tar -xzvf openblas.tar.gz && rm openblas.tar.gz

rm -rf fftw gsl libint openblas && mv fftw* fftw ; mv gsl* gsl ; mv libint* libint ; mv OpenBLAS* openblas

cd fftw     && ./configure --enable-shared --prefix="$PWD/install" && make -j 2 && make                       install && cd ..
cd gsl      && ./configure                 --prefix="$PWD/install" && make -j 2 && make                       install && cd ..
cd openblas &&                                                        make -j 2 && make PREFIX="$PWD/install" install && cd ..

cd libint && cmake -DCMAKE_INSTALL_PREFIX="$PWD/install" -DREQUIRE_CXX_API=OFF . && cmake --build . --parallel 2 && cmake --install . && cd ..

rm -rf external && mkdir external && cp -r fftw/install/* gsl/install/* libint/install/* openblas/install/* external/ && rm -rf fftw gsl libint openblas
