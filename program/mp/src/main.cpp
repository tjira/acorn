#include "linalg.h"
#include "timer.h"
#include <argparse.hpp>
#include <torch/torch.h>

using namespace torch::indexing;

#define TOSTRING(x) #x
#define OCC "ijklmnop"
#define VRT "abcdefgh"

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
}
};

int nthread = 1;

std::vector<std::string> split(const std::string &str, char delim) {
    std::vector<std::string> result; std::stringstream ss(str); std::string line;
    while (getline(ss, line, delim)) {result.push_back (line);} return result;
}

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Moller-Plesset Pertrubation Theory Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-o", "--order").help("-- Order of theory.").default_value(2).scan<'i', int>();
    program.add_argument("-n", "--nthreads").help("-- Number of threads to use.").default_value(nthread).scan<'i', int>();

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now(), mpit = Timer::Now();

    // load the system and integrals in MS basis from disk
    MEASURE("SYSTEM AND ONE-ELECTRON INTEGRALS IN MS BASIS READING: ",
        Matrix Vms = Eigen::LoadMatrix("V_MS.mat");
        Matrix Tms = Eigen::LoadMatrix("T_MS.mat");
        Vector N   = Eigen::LoadMatrix("N.mat"   );
    )
    MEASURE("SYSTEM AND TWO-ELECTRON INTEGRALS IN MS BASIS READING: ",
        torch::Tensor Jms = torch::from_blob(Eigen::LoadTensor("J_MS.mat").data(), {Vms.rows(), Vms.rows(), Vms.rows(), Vms.rows()}, torch::TensorOptions().dtype(torch::kDouble)).clone();
        torch::Tensor Ems = torch::from_blob(Eigen::LoadMatrix("E_MS.mat").data(), {Vms.rows()                                    }, torch::TensorOptions().dtype(torch::kDouble)).clone();
    )

    // extract the number of occupied and virtual orbitals and define the energy
    int nocc = N(0); int nvirt = Vms.rows() / 2 - nocc; int nos = 2 * nocc, nvs = 2 * nvirt; double E = 0; std::vector<double> Empn; nthread = program.get<int>("-n");

    // initialize the antisymmetrized Coulomb integrals in Physicist's notation and the Hamiltonian matrix in MS basis
    torch::Tensor Jmsa = (Jms - Jms.permute({0, 3, 2, 1})).permute({0, 2, 1, 3}); Matrix Hms = Tms + Vms;

    // calculate the HF energy
    for (int i = 0; i < 2 * nocc; i++) {
        E += Hms(i, i); for (int j = 0; j < 2 * nocc; j++) E += 0.5 * Jmsa.index({i, j, i, j}).item().toDouble();
    }

    // print new line
    std::cout << std::endl;

    // loop over the orders of MP perturbation theory
    for (int i = 2; i <= program.get<int>("-o"); i++) {

        // start the timer and define contractions with energy
        mpit = Timer::Now(); std::vector<std::string> contributions = split(mbptf.at(i - 2), ' '); std::vector<double> Empi(contributions.size());

        // print the required number of contributions
        std::printf("MP%02d ENERGY CALCULATION\n%13s %17s %12s\n", i, "CONTR", "VALUE", "TIME");

        // loop over the contractions
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(nthread)
        #endif
        for (int j = 0; j < (int)contributions.size(); j++) {

            // start the timer, split the contribution and replace the dashes with commas
            tp = Timer::Now(); std::vector<std::string> contribution = split(contributions.at(j), ';'); std::replace(contribution.at(0).begin(), contribution.at(0).end(), '-', ',');
            
            // define the tensor slices with axes, array of orbital energy slices and vector of ones for reshaping
            std::vector<torch::Tensor> views; std::vector<Slice> axes(4); std::vector<torch::Tensor> epsts; std::vector<int64_t> ones = {-1};

            // extract the contraction indices and split them
            std::vector<std::string> contrinds = split(contribution.at(0), ',');

            // add the coulomb contributions
            for (int k = 0; k < i; k++) {
                for (int l = 0; l < 4; l++) {
                    if (std::string(OCC).find(contrinds.at(k).at(l)) != std::string::npos) axes.at(l) = Slice(None, nos);
                    if (std::string(VRT).find(contrinds.at(k).at(l)) != std::string::npos) axes.at(l) = Slice(nos, None);
                } views.push_back(Jmsa.index({axes.at(0), axes.at(1), axes.at(2), axes.at(3)}));
            }

            // add the energy denominators
            for (int k = i; k < contrinds.size(); k++, epsts.clear(), ones = {-1}) {
                for (int l = 0; l < contrinds.at(k).size(); l++) {
                    if (std::string(OCC).find(contrinds.at(k).at(l)) != std::string::npos) {
                        epsts.push_back(+1*Ems.index({Slice(None, nos)}).reshape(at::IntArrayRef(ones))), ones.push_back(1);
                    }
                }
                for (int l = 0; l < contrinds.at(k).size(); l++) {
                    if (std::string(VRT).find(contrinds.at(k).at(l)) != std::string::npos) {
                        epsts.push_back(-1*Ems.index({Slice(nos, None)}).reshape(at::IntArrayRef(ones))), ones.push_back(1);
                    }
                }
                views.push_back(1 / (std::accumulate(epsts.begin() + 1, epsts.end(), epsts.at(0))));
            }

            // add the contribution to the energy
            double Empii = std::stod(contribution.at(1)) * torch::einsum(contribution.at(0), views).item().toDouble(); Empi.at(j) = Empii;

            // print the contribution
            std::printf("%06d/%06d %17.14f %s\n", j + 1, (int)contributions.size(), Empii, Timer::Format(Timer::Elapsed(tp)).c_str());
        }

        // print the MPN energy
        Empn.push_back(std::accumulate(Empi.begin(), Empi.end(), 0.0)); std::printf("MP%02d CORRELATION ENERGY: %.16f\n", i, Empn.back());

        // print the MPN timing
        std::printf("MP%02d CORRELATION ENERGY TIMING: %s\n", i, Timer::Format(Timer::Elapsed(mpit)).c_str());

        // print new line
        if (i < program.get<int>("-o")) std::cout << std::endl;
    }

    // sum the MPN energies
    E += std::accumulate(Empn.begin(), Empn.end(), 0.0);

    // print the final energy
    std::printf("\nFINAL SINGLE POINT ENERGY: %.13f\n\nTOTAL TIME: %s\n", E + N(1), Timer::Format(Timer::Elapsed(start)).c_str());
}
