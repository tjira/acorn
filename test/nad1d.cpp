#include "modelsolver.h"

constexpr double result = 0.9, precision = 1e-6;

int test_nad1d(int, char**) {
    // create the model system
    ModelSystem system(2000, {{"0.01*tanh(0.6*x)", "0.001*exp(-x^2)"}, {"0.001*exp(-x^2)", "-0.01*tanh(0.6*x)"}}, {"x"}, {-24, 24}, 4096);

    // create the classical solver options with some parameters
    ModelSolver::OptionsDynamics opt; opt.iters = 350, opt.momentum = {10.95}, opt.position = {-10}, opt.state = 1; opt.seed = 1, opt.step = 10;
    opt.trajs = 10, opt.gradient = {{"0.006/cosh(0.6*x)^2", "-0.002*x*exp(-x^2)"}, {"-0.002*x*exp(-x^2)", "-0.006/cosh(0.6*x)^2"}};
    opt.adiabatic = false, opt.savetraj = false;

    // create the solver and perform the dynamics
    Result res = ModelSolver(opt).run(system, {}, false);

    // print the population of the S0 state and the difference from the reference result
    std::printf("S0 POPULATION: %.14f, DIFFERENCE: %.3e\n", res.msv.pops(0), std::abs(res.msv.pops(0) - result));

    // return the status of the test
    return std::abs(res.msv.pops(0) - result) < precision ? 0 : 1;
}
