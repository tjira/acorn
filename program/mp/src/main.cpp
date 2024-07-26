#include "linalg.h"
#include "timer.h"
#include <argparse.hpp>

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

double loop(const std::vector<std::pair<Tensor<4>, std::vector<int>>>& Jcs, const std::pair<Vector, std::vector<std::vector<int>>>& Ecs, const std::vector<int> dims, std::vector<int> index = {}, size_t i = 0) {
    size_t rank = dims.size(), nos = dims.front(), nvs = dims.back(); double E = 1, D = 0; if (i == 0) index.resize(rank);

    if (i != rank) {
        for (int j = 0; i < rank && j < dims.at(i); j++) {
            index.at(i) = j, E += loop(Jcs, Ecs, dims, index, i + 1);
        }
    } else {
        for (const std::pair<Tensor<4>, std::vector<int>>& contr : Jcs) {
            E *= contr.first(index.at(contr.second.at(0)), index.at(contr.second.at(1)), index.at(contr.second.at(2)), index.at(contr.second.at(3)));
        }
        for (size_t j = 0; j < Ecs.second.size(); j++) {
            for (size_t k = 0, l = Ecs.second.at(j).size(); k < l / 2; k++) D -= Ecs.first(nos + index.at(Ecs.second.at(j).at(k)));
            for (size_t k = 0, l = Ecs.second.at(j).size(); k < l / 2; k++) D += Ecs.first(index.at(Ecs.second.at(j).at(l/2 + k)));
            E /= D, D = 0;
        }
        return E;
    }

    return E - 1;
}

std::vector<std::string> split(const std::string &str, char delim) {
    std::vector<std::string> result; std::stringstream ss(str); std::string line;
    while (getline(ss, line, delim)) {result.push_back (line);} return result;
}

std::tuple<std::vector<std::pair<Tensor<4>, std::vector<int>>>, std::vector<std::vector<int>>, std::vector<int>> parse(const Tensor<4>& Jmsa, const std::string& contr, int order, int nos, int nvs) {
    // resulting slices, contraction axes and dimension
    std::vector<std::pair<Tensor<4>, std::vector<int>>> Jcs;
    std::vector<std::vector<int>> Ecs; std::vector<int> dim;

    // unique indices, index to axis map, and index counter
    std::set<char> indices; std::map<char, int> i2a; int i = 0; Eigen::array<Eigen::Index, 4> si, sr;

    std::vector<std::string> contrs = split(contr, '-');

    // fill the unique indices
    for (const char& c : contr) if (c != '-') indices.insert(c);

    // fill the index to axis map
    for (const char& index : indices) if (std::string(OCC).find(index) != std::string::npos) i2a[index] = i++, dim.push_back(nos);
    for (const char& index : indices) if (std::string(VRT).find(index) != std::string::npos) i2a[index] = i++, dim.push_back(nvs);

    // fill the slices and contraction axes
    for (int j = 0; j < (int)contrs.size(); j++) {
        std::vector<int> axis; for (const char& c : contrs.at(j)) axis.push_back(i2a[c]);
        for (int k = 0; j < order && k < 4; k++) {
            if (std::string(OCC).find(contrs.at(j).at(k)) != std::string::npos) {si[k] = 0, sr[k] = nos;} else {si[k] = nos, sr[k] = nvs;}
        } if (j < order) Jcs.push_back({Jmsa.slice(si, sr), axis}); else Ecs.push_back(axis);
    }

    return {Jcs, Ecs, dim};
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
    MEASURE("SYSTEM AND INTEGRALS IN MS BASIS READING: ",
        Matrix    Vms = Eigen::LoadMatrix("V_MS.mat");
        Matrix    Tms = Eigen::LoadMatrix("T_MS.mat");
        Tensor<4> Jms = Eigen::LoadTensor("J_MS.mat");
        Matrix    Ems = Eigen::LoadMatrix("E_MS.mat");
        Vector    N   = Eigen::LoadMatrix("N.mat"   );
    )

    // extract the number of occupied and virtual orbitals and define the energy
    int nocc = N(0); int nvirt = Jms.dimension(0) / 2 - nocc; int nos = 2 * nocc, nvs = 2 * nvirt; double E = 0; std::vector<double> Empn; nthread = program.get<int>("-n");

    // initialize the antisymmetrized Coulomb integrals in Physicist's notation and the Hamiltonian matrix in MS basis
    Tensor<4> Jmsa = (Jms - Jms.shuffle(Eigen::array<int, 4>{0, 3, 2, 1})).shuffle(Eigen::array<int, 4>{0, 2, 1, 3}); Matrix Hms = Tms + Vms;

    // calculate the HF energy
    for (int i = 0; i < 2 * nocc; i++) {
        E += Hms(i, i); for (int j = 0; j < 2 * nocc; j++) E += 0.5 * Jmsa(i, j, i, j);
    }

    // print new line
    std::cout << std::endl;

    // loop over the orders of MP perturbation theory
    for (int o = 2; o <= program.get<int>("-o"); o++) {

        // start the timer and define contractions with energy
        mpit = Timer::Now(); std::vector<std::string> contrs = split(mbptf.at(o - 2), ' '); std::vector<double> Empi(contrs.size());

        // print the required number of contributions
        std::printf("MP%02d ENERGY CALCULATION\n%13s %17s %12s\n", o, "CONTR", "VALUE", "TIME");

        // loop over the contractions
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(nthread)
        #endif
        for (int i = 0; i < (int)contrs.size(); i++) {

            // start the timer
            tp = Timer::Now();

            // parse the contraction to get slices, axes and dimensions
            auto [Jcs, Ecs, dim] = parse(Jmsa, split(contrs.at(i), ';').at(0), o, nos, nvs);

            // add the contribution to the energy
            double Empii = std::stof(split(contrs.at(i), ';').at(1)) * loop(Jcs, {Ems, Ecs}, dim); Empi.at(i) = Empii;

            // print the contribution
            std::printf("%06d/%06d %17.14f %s\n", i + 1, (int)contrs.size(), Empii, Timer::Format(Timer::Elapsed(tp)).c_str());
        }

        // print the MPN energy
        Empn.push_back(std::accumulate(Empi.begin(), Empi.end(), 0.0)); std::printf("MP%02d CORRELATION ENERGY: %.16f\n", o, Empn.back());

        // print the MPN timing
        std::printf("MP%02d CORRELATION ENERGY TIMING: %s\n", o, Timer::Format(Timer::Elapsed(mpit)).c_str());

        // print new line
        if (o < program.get<int>("-o")) std::cout << std::endl;
    }

    // sum the MPN energies
    E += std::accumulate(Empn.begin(), Empn.end(), 0.0);

    // print the final energy
    std::printf("\nFINAL SINGLE POINT ENERGY: %.13f\n\nTOTAL TIME: %s\n", E + N(1), Timer::Format(Timer::Elapsed(start)).c_str());
}
