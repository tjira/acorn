#pragma once

#include "propagator.h"
#include "result.h"
#include "timer.h"

class ModelSolver {
public:
    struct OptionsAdiabatic {OptionsAdiabatic(){};
        // options structures
        struct Optimize {
            double step=1; int iters=1000;
        } optimize={};
        struct Spectrum {
            bool normalize=false, zpesub=false; int zeropad=0; std::string potential="", window="";
        } spectrum={};

        // variables
        int nstate=1, iters=1000; double step=0.1;
        bool real=false, savewfn=false;
        std::string guess;
    };
    struct OptionsNonadiabatic {OptionsNonadiabatic(){};
        // variables
        int iters=1000; std::vector<std::string> guess; std::string cap="0";
        bool savewfn=false, adiabatic=true; double step=0.1, momentum=0;
    };
    struct OptionsDynamics {
        int iters=100, seed=1, state=0, trajs=1; double step=1.0;
        std::vector<double> position, momentum;
        std::vector<std::string> gradient;
    } dynamics={};

public:
    // constructors and destructors
    ModelSolver(const OptionsNonadiabatic& optn) : optn(optn) {};
    ModelSolver(const OptionsAdiabatic& opta) : opta(opta) {};
    ModelSolver(const OptionsDynamics& optd) : optd(optd) {};

    // methods and solvers
    Result run(const ModelSystem& system, Result res = {}, bool print = true) const;

private:
    Result runnad(const ModelSystem& system, Result res = {}, bool print = true) const;
    Result runad(const ModelSystem& system, Result res = {}, bool print = true) const;
    Result runcd(const ModelSystem& system, Result res = {}, bool print = true) const;

private:
    OptionsAdiabatic opta; OptionsNonadiabatic optn; OptionsDynamics optd;
};
