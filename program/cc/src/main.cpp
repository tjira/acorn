#include "timer.h"
#include "cc.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Moller-Plesset Pertrubation Theory Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now(), tp;

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

    // load the system and integrals in MS basis from disk
    MEASURE("SYSTEM AND INTEGRALS IN MS BASIS READING: ",
        torch::Tensor Jms = torch::ReadTensor("J_MS.mat").squeeze();
        torch::Tensor Vms = torch::ReadTensor("V_MS.mat").squeeze();
        torch::Tensor Tms = torch::ReadTensor("T_MS.mat").squeeze();
        torch::Tensor Ems = torch::ReadTensor("E_MS.mat").squeeze();
        torch::Tensor Fms = torch::ReadTensor("F_MS.mat").squeeze();
        torch::Tensor N   = torch::ReadTensor("N.mat"   ).squeeze();
    )

    // extract the number of occupied and virtual orbitals
    int nocc = N.index({0}).item<int>(); int nvirt = Ems.sizes().at(0) / 2 - nocc; int nos = 2 * nocc, nvs = 2 * nvirt; double E = 0, Ecc = 0, Eccp = 0;

    // extract the command line parameters and define "o" and "v" index ranges
    std::vector<int> exc = program.get<std::vector<int>>("-e"), pts = program.get<std::vector<int>>("-p"); bool linear = program.get<bool>("--linear");
    int iters = program.get<int>("-i"); double thresh = program.get<double>("-t"); auto o = Slice(None, nos), v = Slice(nos, None);

    // initialize the antisymmetrized Coulomb integrals in Physicists' notation and the Hamiltonian matrix in MS basis
    torch::Tensor Jmsa = (Jms - Jms.permute({0, 3, 2, 1})).permute({0, 2, 1, 3}); torch::Tensor Hms = Vms + Tms;

    // create the tripla excitation orbital energy tensor
    torch::Tensor Emst = 1 / (Ems.index({o}).reshape({-1}) + Ems.index({o}).reshape({-1, 1}) + Ems.index({o}).reshape({-1, 1, 1}) - Ems.index({v}).reshape({-1, 1, 1, 1}) - Ems.index({v}).reshape({-1, 1, 1, 1, 1}) - Ems.index({v}).reshape({-1, 1, 1, 1, 1, 1}));

    // create the double excitation orbital energy tensor
    torch::Tensor Emsd = 1 / (Ems.index({o}).reshape({-1}) + Ems.index({o}).reshape({-1, 1}) - Ems.index({v}).reshape({-1, 1, 1}) - Ems.index({v}).reshape({-1, 1, 1, 1}));

    // create the single excitation orbital energy tensor
    torch::Tensor Emss = 1 / (Ems.index({o}).reshape({-1}) - Ems.index({v}).reshape({-1, 1}));

    // initialize single and double excitation amplitudes
    torch::Tensor T1 = torch::zeros({nvs, nos}, torch::dtype(torch::kDouble)), T2 = Jmsa.index({v, v, o, o}) * Emsd;

    // calculate the HF energy
    for (int i = 0; i < 2 * nocc; i++) {
        E += Hms.index({i, i}).item<double>(); for (int j = 0; j < 2 * nocc; j++) E += 0.5 * Jmsa.index({i, j, i, j}).item<double>();
    }

    // print the required number of contributions
    std::printf("\nCC ENERGY CALCULATION\n%13s %17s %12s\n", "CONTR", "VALUE", "TIME");

    // update amplitudes to self consistency
    for (int i = 0; i < iters; i++) {

        // start the timer
        Timer::Timepoint ait = Timer::Now();

        // update the amplitudes
        if (std::find(exc.begin(), exc.end(), 1) != exc.end() && std::find(exc.begin(), exc.end(), 2) != exc.end() && exc.size() == 2 && !linear) {
            std::tie(T1, T2) = Acorn::CC::CCSD::amplitude(Jmsa, Fms, Emss, Emsd, T1, T2, nos), Eccp = Ecc, Ecc = Acorn::CC::CCSD::energy(Jmsa, Fms, T1, T2, nos);
        } else if (std::find(exc.begin(), exc.end(), 2) != exc.end() && exc.size() == 1 && !linear) {
            T2 = Acorn::CC::CCD::amplitude(Jmsa, Emsd, T2, nos), Eccp = Ecc, Ecc = Acorn::CC::CCD::energy(Jmsa, T2, nos);
        } else if (std::find(exc.begin(), exc.end(), 2) != exc.end() && exc.size() == 1 && linear) {
            T2 = Acorn::CC::LCCD::amplitude(Jmsa, Emsd, T2, nos), Eccp = Ecc, Ecc = Acorn::CC::LCCD::energy(Jmsa, T2, nos);
        } else {
            throw std::runtime_error("NO VALID EXCITATIONS SELECTED FOR COUPLED CLUSTER CALCULATION");
        }

        // print the iteration info
        std::printf("%6d %20.14f %.2e %s\n", i + 1, Ecc, std::abs(Ecc - Eccp), Timer::Format(Timer::Elapsed(ait)).c_str());

        // finish if covergence reached
        if (std::abs(Ecc - Eccp) < thresh) {std::cout << std::endl; break;}
        else if (i == iters - 1) throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN THE CC AMPLITUDE SCF REACHED");
    }

    // calculate the perturbation corrections
    if (std::find(exc.begin(), exc.end(), 1) != exc.end() && std::find(exc.begin(), exc.end(), 2) != exc.end() && exc.size() == 2 && !linear) {
        if (std::find(pts.begin(), pts.end(), 3) != pts.end() && pts.size() == 1) {
            MEASURE("THIRD ORDER EXCITATIONS CORRECTION: ", E += Acorn::CC::CCSD::pertrubationTriple(Jmsa, Emst, T1, T2, nos)) std::cout << std::endl;
        } else if (pts.size() > 1) {
            throw std::runtime_error("NO VALID PERTURBATIONS SELECTED FOR COUPLED CLUSTER CALCULATION");
        }
    } else if (pts.size() > 0) {
        throw std::runtime_error("NO VALID EXCITATIONS SELECTED FOR COUPLED CLUSTER CALCULATION");
    }

    // print the final energy
    std::printf("FINAL SINGLE POINT ENERGY: %.13f\n\nTOTAL TIME: %s\n", E + Ecc + N.index({1}).item<double>(), Timer::Format(Timer::Elapsed(start)).c_str());
}
