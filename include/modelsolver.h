#pragma once

#include "modelsystem.h"
#include "result.h"
#include "timer.h"

class ModelSolver {
public:
    struct OptionsAdiabatic {OptionsAdiabatic(){};
        // options structures
        struct Dynamics {
            int iters=100, state=0; double step=1.0;
            std::vector<double> position, velocity;
            std::vector<std::string> gradient;
        } dynamics={};
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
        // dynamics structure
        struct Dynamics {
            int iters=100, state=0; double step=1.0;
            std::vector<double> position, velocity;
            std::vector<std::string> gradient;
        } dynamics={};

        // variables
        int iters=1000; std::vector<std::string> guess; std::string cap="0";
        bool savewfn=false, adiabatic=true; double step=0.1, momentum=0;
    };
public:
    // constructors and destructors
    ModelSolver(const OptionsNonadiabatic& optn) : optn(optn) {};
    ModelSolver(const OptionsAdiabatic& opta) : opta(opta) {};

    // methods and solvers
    Result run(const ModelSystem& system, Result res = {}, bool print = true);

private:
    template <typename T> Result runcd(const ModelSystem& system, const T& optdyn, Result res = {}, bool print = true);
    Result runnad(const ModelSystem& system, Result res = {}, bool print = true);
    Result runad(const ModelSystem& system, Result res = {}, bool print = true);

private:
    OptionsAdiabatic opta; OptionsNonadiabatic optn;
};
