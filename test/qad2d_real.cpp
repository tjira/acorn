#include "modelsolver.h"

constexpr double result = 1.25157958326112, precision = 1e-6;

int test_qad2d_real(int, char**) {
    // create the model system
    ModelSystem system(1, {{"0.5*(x^2+y^2)"}}, {"x", "y"}, {-16, 16}, 256);

    // create the adiabatic solver options with guess function
    ModelSolver::OptionsAdiabatic opt; opt.guess = "exp(-(x^2+y^2))", opt.iters = 200;
    opt.nstate = 1, opt.real = true, opt.savewfn = false, opt.step = 0.1;

    // create the solver and perforn the dynamics
    Result res = ModelSolver(opt).run(system, {}, false);

    // print the total energy and the difference from the reference result
    std::printf("ENERGY: %.14f, DIFFERENCE: %.3e\n", res.msv.energy(0), std::abs(res.msv.energy(0) - result));

    // return the status of the test
    return std::abs(res.msv.energy(0) - result) < precision ? 0 : 1;
}
