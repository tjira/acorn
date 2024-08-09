#pragma once

#include <string>
#include <vector>

namespace Acorn {
    namespace CDYN {
        struct Options {
            std::vector<std::string> potential; double mass, momentum, position, step; int excstate, iters, log, nstate, seed, trajs; bool adiabatic, savetraj;
        };
    }
}
