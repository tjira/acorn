#include "acorn.h"

int main(int argc, char** argv) {
    // get executable path for the executable to run from anywhere
    ip = std::filesystem::weakly_canonical(argv[0]).parent_path();

    // create the model system
    ModelSystem system(1, {{"0.5*x^2"}}, {-16, 16}, 1024);

    // create the adiabatic solver options with mandatory values
    ModelSolver::OptionsAdiabatic opta; opta.guess = "exp(-x^2)", opta.savewfn = true;

    // create the solver and perforn the dynamics
    Result res = ModelSolver(opta).run(system, {}, false);

    // print the total energy
    std::printf("%.8f\n", res.msv.opten(0));
}
