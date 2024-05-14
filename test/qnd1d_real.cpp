#include "modelsolver.h"

constexpr double result = 0.02022577692910, precision = 1e-6;

int test_qnd1d_real(int, char**) {
    // create the model system
    ModelSystem system(2000, {{"0.01*tanh(0.6*x)", "0.001*exp(-x^2)"}, {"0.001*exp(-x^2)", "-0.01*tanh(0.6*x)"}}, {"x"}, {-24, 24}, 4096);

    // create the nonadiabatic solver options with some parameters
    ModelSolver::OptionsNonadiabatic opt; opt.adiabatic = false, opt.guess = {"exp(-(x+10)^2)"},
    opt.iters = 350, opt.momentum = 10.95, opt.savewfn = false, opt.step = 10, opt.cap = "0";

    // create the solver and perform the dynamics
    Result res = ModelSolver(opt).run(system, {}, false);

    // print the total energy and the difference from the reference result
    std::printf("ENERGY: %.14f, DIFFERENCE: %.3e\n", res.msv.energy(0), std::abs(res.msv.energy(0) - result));

    // return the status of the test
    return std::abs(res.msv.energy(0) - result) < precision ? 0 : 1;
}
