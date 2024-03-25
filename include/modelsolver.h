#pragma once

#include "modelsystem.h"
#include "result.h"

class ModelSolver {
public:
    struct OptionsAdiabatic {OptionsAdiabatic(){};
        // options structures
        struct Dynamics {
            struct Berendsen {double tau=1, temp=0, timeout=10;} berendsen={};
            int iters=100; double step=1.0;
        } dynamics={};
        struct Spectrum {
            bool zpesub=false; std::string potential="", window="";
        } spectrum={};

        // variables
        bool real=false, optimize=false, savewfn=false;
        int nstate=1, iters=1000; double step=0.1;
        std::string guess;
    };
    struct OptionsNonadiabatic {OptionsNonadiabatic(){};
        // dynamics structure
        struct Dynamics {
            struct Berendsen {double tau=1, temp=0, timeout=10;} berendsen={};
            int iters=100; double step=1.0;
        } dynamics={};

        // variables
        bool savewfn=false; double step=0.1; int iters=1000;
        std::vector<std::string> guess;
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
