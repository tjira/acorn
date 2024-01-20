#pragma once

#include "modelsystem.h"

class ModelSolver {
public:
    struct Options {Options(){};
        // structures
        struct Dynamics {int iters=100; double step=1.0; std::string output="trajectory.xyz";} dynamics={};

        // variables
        int nstate=3, iters=1000; bool real=false, optimize=false; double step=0.1, thresh=1e-8;
    };
    struct Result {
        std::vector<std::vector<Vector<std::complex<double>>>> states;
        Vector<> energy, r, u;
    };
public:
    // constructors and destructors
    ModelSolver(const Options& opt) : opt(opt) {};

    // methods and solvers
    Result run(const ModelSystem& system, Result res = {}, bool print = true);

private:
    Options opt;
};
