#include "mp.h"
#include "timer.h"
#include <argparse.hpp>

#define TOSTRING(x) #x
std::vector<std::string> mbptf {
std::string {
#include "mbpt2.txt"
},
std::string {
#include "mbpt3.txt"
},
std::string {
#include "mbpt4.txt"
},
std::string {
#include "mbpt5.txt"
},
std::string {
#include "mbpt6.txt"
},
};

int nthread = 1;

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Moller-Plesset Perturbation Theory Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now(), tp;

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-o", "--order").help("-- Order of the theory.").default_value(2).scan<'i', int>();
    program.add_argument("-n", "--nthreads").help("-- Number of threads to use.").default_value(nthread).scan<'i', int>();

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} std::vector<Timer::Timepoint> timers(2);

    // set the number of threads
    nthread = program.get<int>("-n");

    // extract the command line parameters
    Acorn::MBPT::Options opt = {program.get<int>("-o")};

    // load the system and integrals in MS basis from disk
    MEASURE("SYSTEM AND INTEGRALS IN MS BASIS READING: ",
        torch::Tensor Jms = torch::ReadTensor("J_MS.mat").squeeze();
        torch::Tensor Vms = torch::ReadTensor("V_MS.mat").squeeze();
        torch::Tensor Tms = torch::ReadTensor("T_MS.mat").squeeze();
        torch::Tensor Ems = torch::ReadTensor("E_MS.mat").squeeze();
        torch::Tensor N   = torch::ReadTensor("N.mat"   ).squeeze();
    )

    // extract the number of occupied and virtual spinorbitals
    int nos = 2 * N.index({0}).item<int>(); int nvs = Ems.sizes().at(0) - nos;

    // define the energy containers
    std::vector<double> Empn; double E = 0;

    // initialize the antisymmetrized Coulomb integrals in Physicists' notation and the Hamiltonian matrix in MS basis
    torch::Tensor Jmsa = (Jms - Jms.permute({0, 3, 2, 1})).permute({0, 2, 1, 3}); torch::Tensor Hms = Vms + Tms;

    // calculate the HF energy
    for (int i = 0; i < nos; i++) {
        E += Hms.index({i, i}).item<double>(); for (int j = 0; j < nos; j++) E += 0.5 * Jmsa.index({i, j, i, j}).item<double>();
    }

    // print new line
    std::cout << std::endl;

    // loop over the orders of MP perturbation theory
    for (int i = 2; i <= program.get<int>("-o"); i++) {

        // start the timer
        timers.at(0) = Timer::Now();

        // define contributions with energy
        std::vector<std::string> contributions = SPLIT(mbptf.at(i - 2), ' '); std::vector<double> Empi(contributions.size());

        // print the required number of contributions
        std::printf("MP%02d ENERGY CALCULATION\n%13s %17s %12s\n", i, "CONTR", "VALUE", "TIME");

        // loop over the contractions
        #pragma omp parallel for num_threads(nthread)
        for (size_t j = 0; j < contributions.size(); j++) {

            // start the timer
            timers.at(1) = Timer::Now();

            // add the contribution to the energy
            double Empii = Acorn::MBPT::evaluate(contributions.at(j), Jmsa, Ems, nos, i); Empi.at(j) = Empii;

            // print the contribution
            std::printf("%06d/%06d %17.14f %s\n", (int)j + 1, (int)contributions.size(), Empii, Timer::Format(Timer::Elapsed(timers.at(1))).c_str());
        }

        // print the MPN energy
        Empn.push_back(std::accumulate(Empi.begin(), Empi.end(), 0.0)); std::printf("MP%02d CORRELATION ENERGY: %.16f\n", i, Empn.back());

        // print the MPN timing
        std::printf("MP%02d CORRELATION ENERGY TIMING: %s\n", i, Timer::Format(Timer::Elapsed(timers.at(0))).c_str());

        // print new line
        if (i < program.get<int>("-o")) std::cout << std::endl;
    }

    // sum the MPN energies
    E += std::accumulate(Empn.begin(), Empn.end(), 0.0);

    // print the final energy
    std::printf("\nFINAL SINGLE POINT ENERGY: %.13f\n\nTOTAL TIME: %s\n", E + N.index({1}).item<double>(), Timer::Format(Timer::Elapsed(start)).c_str());
}
