#!/bin/bash

cat << EndOfMessage
#pragma once
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/FFT>
#include <bits/stdc++.h>
EndOfMessage

tail -n +5 include/constant.h

cat \
    include/eigen.h \
    include/table.h \
    include/system.h \
    include/integral.h \
    include/population.h \
    include/timer.h \
    include/result.h \
    include/method.h  \
    include/restrictedmethod.h  \
    include/unrestrictedmethod.h  \
    include/restrictedhartreefock.h  \
    include/unrestrictedhartreefock.h  \
    include/numpy.h \
    include/transform.h \
    include/determinant.h \
    include/restrictedconfigurationinteraction.h  \
    include/restrictedmollerplesset.h  \
    include/expression.h \
    include/modelsystem.h \
    include/modelsolver.h \
    include/orca.h \
| sed '/#pragma once/d ; /#include/d ; /\/\//d ; /^$/d'
