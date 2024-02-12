#pragma once

#include "modelsystem.h"
#include "result.h"

class ModelSolver {
public:
    struct OptionsAdiabatic {OptionsAdiabatic(){};
        // dynamics structure
        struct Dynamics {
            struct Berendsen {double tau=1, temp=0, timeout=10;} berendsen={};
            int iters=100; double step=1.0; std::string folder;
        } dynamics={};

        // variables
        int nstate=3, iters=1000; bool real=false, optimize=false;
        std::string folder, guess; double step=0.1, thresh=1e-8;
    };
    struct OptionsNonadiabatic {OptionsNonadiabatic(){};
        // dynamics structure
        struct Dynamics {
            struct Berendsen {double tau=1, temp=0, timeout=10;} berendsen={};
            int iters=100; double step=1.0; std::string folder;
        } dynamics={};

        // variables
        std::vector<std::string> guess; std::string folder; double step=0.1; int iters=1000;
    };
public:
    // constructors and destructors
    ModelSolver(const OptionsNonadiabatic& optn) : optn(optn) {};
    ModelSolver(const OptionsAdiabatic& opta) : opta(opta) {};

    // methods and solvers
    Result run(const ModelSystem& system, Result res = {}, bool print = true);

private:
    Result runnad(const ModelSystem& system, Result res = {}, bool print = true);
    Result runad(const ModelSystem& system, Result res = {}, bool print = true);

private:
    OptionsAdiabatic opta; OptionsNonadiabatic optn;
};
