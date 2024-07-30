#include "tensor.h"
#include "timer.h"
#include <argparse.hpp>

using namespace torch::indexing;

#define OCC "abcdefghijklmnopqrstuvwxyz"
#define VRT "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

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

std::vector<std::string> split(const std::string &str, char delim) {
    std::vector<std::string> result; std::stringstream ss(str); std::string line;
    while (getline(ss, line, delim)) {result.push_back (line);} return result;
}

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Moller-Plesset Pertrubation Theory Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now(), tp;

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-o", "--order").help("-- Order of theory.").default_value(2).scan<'i', int>();
    program.add_argument("-n", "--nthreads").help("-- Number of threads to use.").default_value(nthread).scan<'i', int>();

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // set the number of threads
    nthread = program.get<int>("-n");

    // load the system and integrals in MS basis from disk
    MEASURE("SYSTEM AND INTEGRALS IN MS BASIS READING: ",
        torch::Tensor Jms = torch::ReadTensor("J_MS.mat").squeeze();
        torch::Tensor Vms = torch::ReadTensor("V_MS.mat").squeeze();
        torch::Tensor Tms = torch::ReadTensor("T_MS.mat").squeeze();
        torch::Tensor Ems = torch::ReadTensor("E_MS.mat").squeeze();
        torch::Tensor N   = torch::ReadTensor("N.mat"   ).squeeze();
    )

    // extract the number of occupied and virtual orbitals
    int nocc = N.index({0}).item<int>(); int nvirt = Ems.sizes().at(0) / 2 - nocc; int nos = 2 * nocc, nvs = 2 * nvirt;

    // define the energy containers
    std::vector<double> Empn; double E = 0;

    // initialize the antisymmetrized Coulomb integrals in Physicist's notation and the Hamiltonian matrix in MS basis
    torch::Tensor Jmsa = (Jms - Jms.permute({0, 3, 2, 1})).permute({0, 2, 1, 3}); torch::Tensor Hms = Vms + Tms;

    // calculate the HF energy
    for (int i = 0; i < 2 * nocc; i++) {
        E += Hms.index({i, i}).item<double>(); for (int j = 0; j < 2 * nocc; j++) E += 0.5 * Jmsa.index({i, j, i, j}).item<double>();
    }

    // print new line
    std::cout << std::endl;

    // loop over the orders of MP perturbation theory
    for (int i = 2; i <= program.get<int>("-o"); i++) {

        // start the timer
        Timer::Timepoint mpit = Timer::Now();

        // define contributions with energy
        std::vector<std::string> contributions = split(mbptf.at(i - 2), ' '); std::vector<double> Empi(contributions.size());

        // print the required number of contributions
        std::printf("MP%02d ENERGY CALCULATION\n%13s %17s %12s\n", i, "CONTR", "VALUE", "TIME");

        // loop over the contractions
        #pragma omp parallel for num_threads(nthread)
        for (int j = 0; j < (int)contributions.size(); j++) {

            // start the timer
            Timer::Timepoint mpiit = Timer::Now();

            // define the contribution and split it
            std::vector<std::string> contribution = split(contributions.at(j), ';');

            // replace the plus sign with a comma
            std::replace(contribution.at(0).begin(), contribution.at(0).end(), '+', ',');
            std::replace(contribution.at(1).begin(), contribution.at(1).end(), '+', ',');
            
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

            // define the contraction path
            std::vector<std::tuple<std::string, std::array<int, 2>>> path;

            // split the nodes and add them to the path
            for (const std::string& node: split(contribution.at(1), ':')) {
                path.push_back({split(node, '/').at(1), {std::stoi(split(split(node, '/').at(0), '-').at(0)), std::stoi(split(split(node, '/').at(0), '-').at(1))}});
            }

            // contract over every node
            for (const auto& node : path) {

                // contract the nodes
                views.push_back(torch::einsum(std::get<0>(node), {views.at(std::get<1>(node).at(0)), views.at(std::get<1>(node).at(1))}));

                // erase the contracted nodes
                views.erase(views.begin() + std::get<1>(node).at(std::get<1>(node).at(0) < std::get<1>(node).at(1)));
                views.erase(views.begin() + std::get<1>(node).at(std::get<1>(node).at(0) > std::get<1>(node).at(1)));
            }

            // add the contribution to the energy
            double Empii = std::stod(contribution.at(2)) * views.at(0).item<double>(); Empi.at(j) = Empii;

            // print the contribution
            std::printf("%06d/%06d %17.14f %s\n", j + 1, (int)contributions.size(), Empii, Timer::Format(Timer::Elapsed(mpiit)).c_str());
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
    std::printf("\nFINAL SINGLE POINT ENERGY: %.13f\n\nTOTAL TIME: %s\n", E + N.index({1}).item<double>(), Timer::Format(Timer::Elapsed(start)).c_str());
}
