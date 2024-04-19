#include "../bin/acorn.h"

int main(int argc, char** argv) {
    // create the model system
    ModelSystem system(1, {{argv[1]}}, {"x"}, {-16, 16}, 1024);

    // create the adiabatic solver options with guess function
    ModelSolver::OptionsAdiabatic opta; opta.guess = "exp(-x^2)";

    // create the solver and perforn the dynamics
    Result res = ModelSolver(opta).run(system, {}, false);

    // print the total energy
    std::printf("%.8f\n", res.msv.opten(0));
}
