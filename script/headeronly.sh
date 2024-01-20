#!/bin/bash

cat \
    include/constant.h \
    include/eigen.h \
    include/table.h \
    include/system.h \
    include/integral.h \
    include/timer.h \
    include/result.h \
    include/method.h  \
    include/restrictedmethod.h  \
    include/restrictedhartreefock.h  \
    include/transform.h \
    include/restrictedmollerplesset.h  \
    include/expression.h \
    include/modelsystem.h \
    include/modelsolver.h \
    include/orca.h \
    include/default.h \
| sed '/#pragma once/d ; /#include "/d'
