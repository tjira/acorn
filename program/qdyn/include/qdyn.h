#pragma once

#include "fourier.h"
#include "wavefunction.h"

namespace Acorn {
    namespace QDYN {
        struct Options {
            double factor, mass, momentum, step; int dim, iters, optstates; bool adiabatic, align, imaginary, savewfn;
        };
    }
}
