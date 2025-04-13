#!/bin/bash

wget -O     fftw.tar.gz https://www.fftw.org/fftw-3.3.10.tar.gz
wget -O      gsl.tar.gz https://sunsite.icm.edu.pl/pub/gnu/gsl/gsl-2.8.tar.gz
wget -O openblas.tar.gz https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.29/OpenBLAS-0.3.29.tar.gz

tar -xzvf fftw.tar.gz     && rm     fftw.tar.gz
tar -xzvf gsl.tar.gz      && rm      gsl.tar.gz
tar -xzvf openblas.tar.gz && rm openblas.tar.gz

mv fftw* fftw ; mv gsl* gsl ; mv OpenBLAS* openblas

cd fftw && ./configure --prefix="$PWD/install" && make -j 2 && make install && cd ..
cd gsl  && ./configure --prefix="$PWD/install" && make -j 2 && make install && cd ..

cd openblas && make -j 2 && make PREFIX="$PWD/install" install && cd ..

rm -rf external && mkdir external && cp -r fftw/install/* gsl/install/* openblas/install/* external/ && rm -rf fftw gsl openblas
