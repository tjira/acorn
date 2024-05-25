#include "transform.h"
#include "timer.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Hartree-Fock Integral Transform Engine", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // load the integrals in AO basis from disk
    tp = Timer::Now(); std::cout << "V_AO MATRIX READING: " << std::flush; EigenMatrix<> V = Eigen::LoadMatrix("V_AO.mat"); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "T_AO MATRIX READING: " << std::flush; EigenMatrix<> T = Eigen::LoadMatrix("T_AO.mat"); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "S_AO MATRIX READING: " << std::flush; EigenMatrix<> S = Eigen::LoadMatrix("S_AO.mat"); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "J_AO TENSOR READING: " << std::flush; EigenTensor<> J = Eigen::LoadTensor("J_AO.mat"); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "C_MO MATRIX READING: " << std::flush; EigenMatrix<> C = Eigen::LoadMatrix("C_MO.mat"); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // print a new line
    std::cout << std::endl;

    // one-electron integrals in molecular spatial orbital basis
    tp = Timer::Now(); std::cout << "V_AO MATRIX MO TRANSFORM: " << std::flush; Eigen::Write("V_MO.mat", Transform::SingleSpatial(V, C)); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "T_AO MATRIX MO TRANSFORM: " << std::flush; Eigen::Write("T_MO.mat", Transform::SingleSpatial(T, C)); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "S_AO MATRIX MO TRANSFORM: " << std::flush; Eigen::Write("S_MO.mat", Transform::SingleSpatial(S, C)); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // one-electron integrals in molecular spinorbital basis
    tp = Timer::Now(); std::cout << "V_AO MATRIX MS TRANSFORM: " << std::flush; Eigen::Write("V_MS.mat", Transform::SingleSpin(V, C)); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "T_AO MATRIX MS TRANSFORM: " << std::flush; Eigen::Write("T_MS.mat", Transform::SingleSpin(T, C)); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;
    tp = Timer::Now(); std::cout << "S_AO MATRIX MS TRANSFORM: " << std::flush; Eigen::Write("S_MS.mat", Transform::SingleSpin(S, C)); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // two-electron integrals in molecular spatial orbital basis
    tp = Timer::Now(); std::cout << "J_AO TENSOR MO TRANSFORM: " << std::flush; Eigen::Write("J_MO.mat", Transform::CoulombSpatial(J, C)); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // two-electron integrals in molecular spinorbital basis
    tp = Timer::Now(); std::cout << "J_AO TENSOR MS TRANSFORM: " << std::flush; Eigen::Write("J_MS.mat", Transform::CoulombSpin(J, C)); std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}
