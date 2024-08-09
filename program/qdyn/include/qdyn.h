#pragma once

namespace Acorn {
    namespace QDYN {
        struct Options {
            double factor, mass, momentum, step; int dim, iters, optstates; bool adiabatic, align, imaginary, savewfn;
        };
    }
}
