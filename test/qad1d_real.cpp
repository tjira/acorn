#include "modelsolver.h"

constexpr double result = 1.12493167092487, precision = 1e-6;

int test_qad1d_real(int, char**) {
    // create the model system
    ModelSystem system(1, {{"0.5*x^2"}}, {"x"}, {-16, 16}, 1024);

    // create the adiabatic solver options with guess function
    ModelSolver::OptionsAdiabatic opt; opt.guess = "exp(-(x-1)^2)", opt.real = true;

    // create the solver and perforn the dynamics
    Result res = ModelSolver(opt).run(system, {}, false);

    // print the total energy and the difference from the reference result
    std::printf("ENERGY: %.14f, DIFFERENCE: %.3e\n", res.msv.energy(0), std::abs(res.msv.energy(0) - result));

    // return the status of the test
    return std::abs(res.msv.energy(0) - result) < precision ? 0 : 1;
}
