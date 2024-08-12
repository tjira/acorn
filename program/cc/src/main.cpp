#include "cc.h"
#include <argparse.hpp>

#define FORMAT(T) [&](long ms) {char s[99]; std::sprintf(s, "%02ld:%02ld:%02ld.%03ld", ms / 3600000, ms % 3600000 / 60000, ms % 60000 / 1000, ms % 1000); return std::string(s);}(T)

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Coupled Clusters Program", "1.0", argparse::default_arguments::none);

    // define the timers
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> timers(2, std::chrono::high_resolution_clock().now()); auto tp = timers.at(0);
    auto elapsed = [](auto ms) {return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock().now() - ms).count();};

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-i", "--iterations").help("-- Number of SCF iterations.").default_value(100).scan<'i', int>();
    program.add_argument("-e", "--excitations").help("-- Excitations to include in the CC calculation.").nargs(argparse::nargs_pattern::at_least_one).default_value(std::vector<int>{1, 2}).scan<'i', int>();
    program.add_argument("-p", "--perturbations").help("-- Perturbations to include in the CC calculation.").nargs(argparse::nargs_pattern::any).default_value(std::vector<int>{}).scan<'i', int>();
    program.add_argument("-t", "--threshold").help("-- Convergence threshold.").default_value(1e-8).scan<'g', double>();
    program.add_argument("--linear").help("-- Use linearized coupled clusters").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // extract the command line parameters
    Acorn::CC::Options opt = {program.get<std::vector<int>>("-e"), program.get<std::vector<int>>("-p"), program.get<double>("-t"), program.get<int>("-i"), program.get<bool>("--linear")};

    // start the timer for integral loading
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of integral loading
    std::cout << "SYSTEM AND INTEGRALS IN MS BASIS READING: " << std::flush;

    // load the system and integrals in MS basis from disk
    torch::Tensor Jms = torch::ReadTensor("J_MS.mat").squeeze();
    torch::Tensor Vms = torch::ReadTensor("V_MS.mat").squeeze();
    torch::Tensor Tms = torch::ReadTensor("T_MS.mat").squeeze();
    torch::Tensor Ems = torch::ReadTensor("E_MS.mat").squeeze();
    torch::Tensor Fms = torch::ReadTensor("F_MS.mat").squeeze();
    torch::Tensor N   = torch::ReadTensor("N.mat"   ).squeeze();

    // print the time for integral loading
    std::cout << FORMAT(elapsed(timers.at(1))) << std::endl;

    // extract the number of occupied and virtual spinorbitals
    int nos = 2 * N.index({0}).item<int>(); int nvs = Ems.sizes().at(0) - nos;

    // initialize the energy variables and orbital slices
    double E = 0, Ecc = 0, Eccp = 0; auto o = Slice(None, nos), v = Slice(nos, None);

    // initialize the antisymmetrized Coulomb integrals in Physicists' notation and the Hamiltonian matrix in MS basis
    torch::Tensor Jmsa = (Jms - Jms.permute({0, 3, 2, 1})).permute({0, 2, 1, 3}); torch::Tensor Hms = Vms + Tms;

    // create orbital energy tensors
    torch::Tensor Emst = 1 / (Ems.index({o}).reshape({-1}) + Ems.index({o}).reshape({-1, 1}) + Ems.index({o}).reshape({-1, 1, 1}) - Ems.index({v}).reshape({-1, 1, 1, 1}) - Ems.index({v}).reshape({-1, 1, 1, 1, 1}) - Ems.index({v}).reshape({-1, 1, 1, 1, 1, 1}));
    torch::Tensor Emsd = 1 / (Ems.index({o}).reshape({-1}) + Ems.index({o}).reshape({-1, 1}) - Ems.index({v}).reshape({-1, 1, 1}) - Ems.index({v}).reshape({-1, 1, 1, 1}));
    torch::Tensor Emss = 1 / (Ems.index({o}).reshape({-1}) - Ems.index({v}).reshape({-1, 1}));

    // initialize single and double excitation amplitudes
    torch::Tensor T1 = torch::zeros({nvs, nos}, torch::dtype(torch::kDouble)), T2 = Jmsa.index({v, v, o, o}) * Emsd;

    // calculate the HF energy
    for (int i = 0; i < nos; i++) {
        E += Hms.index({i, i}).item<double>(); for (int j = 0; j < nos; j++) E += 0.5 * Jmsa.index({i, j, i, j}).item<double>();
    }

    // print the required number of contributions
    std::printf("\nCC ENERGY CALCULATION\n%13s %17s %12s\n", "CONTR", "VALUE", "TIME");

    // update amplitudes to self consistency
    for (int i = 0; i < opt.iters; i++) {

        // start the timer
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // update the amplitudes
        if (std::find(opt.exc.begin(), opt.exc.end(), 1) != opt.exc.end() && std::find(opt.exc.begin(), opt.exc.end(), 2) != opt.exc.end() && opt.exc.size() == 2 && !opt.linear) {
            std::tie(T1, T2) = Acorn::CC::CCSD::amplitude(Jmsa, Fms, Emss, Emsd, T1, T2, nos), Eccp = Ecc, Ecc = Acorn::CC::CCSD::energy(Jmsa, Fms, T1, T2, nos);
        } else if (std::find(opt.exc.begin(), opt.exc.end(), 2) != opt.exc.end() && opt.exc.size() == 1 && !opt.linear) {
            T2 = Acorn::CC::CCD::amplitude(Jmsa, Emsd, T2, nos), Eccp = Ecc, Ecc = Acorn::CC::CCD::energy(Jmsa, T2, nos);
        } else if (std::find(opt.exc.begin(), opt.exc.end(), 2) != opt.exc.end() && opt.exc.size() == 1 && opt.linear) {
            T2 = Acorn::CC::LCCD::amplitude(Jmsa, Emsd, T2, nos), Eccp = Ecc, Ecc = Acorn::CC::LCCD::energy(Jmsa, T2, nos);
        } else {
            throw std::runtime_error("NO VALID EXCITATIONS SELECTED FOR COUPLED CLUSTER CALCULATION");
        }

        // print the iteration info
        std::printf("%6d %20.14f %.2e %s\n", i + 1, Ecc, std::abs(Ecc - Eccp), FORMAT(elapsed(timers.at(1))).c_str());

        // finish if covergence reached
        if (std::abs(Ecc - Eccp) < opt.thresh) {std::cout << std::endl; break;}
        else if (i == opt.iters - 1) throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN THE CC AMPLITUDE SCF REACHED");
    }

    // calculate the perturbation corrections
    if (std::find(opt.exc.begin(), opt.exc.end(), 1) != opt.exc.end() && std::find(opt.exc.begin(), opt.exc.end(), 2) != opt.exc.end() && opt.exc.size() == 2 && !opt.linear) {
        if (std::find(opt.pts.begin(), opt.pts.end(), 3) != opt.pts.end() && opt.pts.size() == 1) {
            timers.at(1) = std::chrono::high_resolution_clock().now(); std::cout << "THIRD ORDER EXCITATIONS CORRECTION: " << std::flush; E += Acorn::CC::CCSD::perturbationTriple(Jmsa, Emst, T1, T2, nos); std::cout << FORMAT(elapsed(timers.at(1))) << std::endl << std::endl;
        } else if (opt.pts.size() > 1) {
            throw std::runtime_error("NO VALID PERTURBATIONS SELECTED FOR COUPLED CLUSTER CALCULATION");
        }
    } else if (opt.pts.size() > 0) {
        throw std::runtime_error("NO VALID EXCITATIONS SELECTED FOR COUPLED CLUSTER CALCULATION");
    }

    // print the final energy
    std::printf("FINAL SINGLE POINT ENERGY: %.13f\n\nTOTAL TIME: %s\n", E + Ecc + N.index({1}).item<double>(), FORMAT(elapsed(timers.at(0))).c_str());
}
