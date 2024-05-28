#include "system.h"
#include "timer.h"
#include <argparse.hpp>

inline bool VectorContains(const std::vector<int>& v, const int& e) {return std::find(v.begin(), v.end(), e) != v.end();}

std::tuple<std::array<std::vector<int>, 2>, int> Align(std::array<std::vector<int>, 2> det1, const std::array<std::vector<int>, 2>& det2) {
    // define the temporary determinant and swaps variable
    int swaps = 0;
    
    // align the alpha electrons
    for (size_t i = 0; i < det1.at(0).size(); i++) {
        if (det1.at(0).at(i) != det2.at(0).at(i)) {
            for (size_t j = 0; j < det1.at(0).size(); j++) {
                if (det1.at(0).at(i) == det2.at(0).at(j)) {
                    std::swap(det1.at(0).at(i), det1.at(0).at(j)), swaps++;
                }
            }
        }
    }

    // align the beta electrons
    for (size_t i = 0; i < det1.at(1).size(); i++) {
        if (det1.at(1).at(i) != det2.at(1).at(i)) {
            for (size_t j = 0; j < det1.at(1).size(); j++) {
                if (det1.at(1).at(i) == det2.at(1).at(j)) {
                    std::swap(det1.at(1).at(i), det1.at(1).at(j)), swaps++;
                }
            }
        }
    }

    // return the aligned determinant
    return {det1, swaps};
}

std::vector<std::vector<int>> Combinations(int n, int k) {
    // create the bitmask that will get permuted and the resulting vector
    std::string bitmask(k, 1); bitmask.resize(n, 0); std::vector<std::vector<int>> combs;
 
    // generate the combinations
    do {std::vector<int> comb; comb.reserve(k);
        for (int j = 0; j < n; j++) {
            if (bitmask[j]) comb.push_back(j);
        } combs.push_back(comb);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    // return the result
    return combs;
}

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Configuration Interaction Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-f", "--file").help("-- System file in the .xyz format.").default_value("molecule.xyz");

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // load the integrals in AO basis and system from disk
    MEASURE("SYSTEM AND INTEGRALS IN MS BASIS READING: ",
        EigenMatrix<> Vms = Eigen::LoadMatrix("V_MS.mat");
        EigenMatrix<> Tms = Eigen::LoadMatrix("T_MS.mat");
        EigenTensor<> Jms = Eigen::LoadTensor("J_MS.mat");
        System system(program.get("-f"));
    )

    // define the number of basis functions, Hamiltonian in MS basis and determinant container
    int nbf = Vms.rows() / 2; EigenMatrix<> Hms = Tms + Vms; std::vector<std::array<std::vector<int>, 2>> dets;

    MEASURE("GENERATING ALL POSSIBLE DETERMINANTS:     ",
        for (const std::vector<int>& alpha : Combinations(nbf, system.nocc())) {
            for (const std::vector<int>& beta : Combinations(nbf, system.nocc())) {
                dets.push_back(std::array<std::vector<int>, 2>{alpha, beta});
            }
        }
    ) std::cout << "\nNUMBER OF DETERMINANTS GENERATED: " << dets.size() << std::endl;

    // define the CI Hamiltonian
    EigenMatrix<> Hci = EigenMatrix<>::Zero(dets.size(), dets.size());

    // fill the CI Hamiltonian
    for (int i = 0; i < Hci.rows(); i++) {
        for (int j = 0; j < Hci.cols(); j++) {
            auto det1 = dets.at(i); auto [det2, swaps] = Align(dets.at(j), det1); int diff = 0; double elem = 0;

             // define all needed vectors of spinorbitals
            std::vector<int> so1(2 * det1.at(0).size()), so2(2 * det1.at(0).size()), common, unique;

            // calculate the number of differences
            for (size_t i = 0; i < det1.at(0).size(); i++) if (det1.at(0).at(i) != det2.at(0).at(i)) diff++;
            for (size_t i = 0; i < det1.at(0).size(); i++) if (det1.at(1).at(i) != det2.at(1).at(i)) diff++;

            // fill the spinorbitals of det1 and det2 determinant
            for (size_t i = 0; i < det1.at(0).size(); i++) so2.at(i) = 2 * det2.at(0).at(i), so2.at(det2.at(0).size() + i) = 2 * det2.at(1).at(i) + 1;
            for (size_t i = 0; i < det1.at(0).size(); i++) so1.at(i) = 2 * det1.at(0).at(i), so1.at(det1.at(0).size() + i) = 2 * det1.at(1).at(i) + 1;

            // fill the common spinorbitals vector
            for (size_t i = 0; i < so1.size(); i++) if (so1.at(i) == so2.at(i)) common.push_back(so1.at(i));

            // fill the unique spinorbitals vector
            for (size_t i = 0; i < so1.size(); i++) if (!VectorContains(so1, so2.at(i))) unique.push_back(so2.at(i));
            for (size_t i = 0; i < so1.size(); i++) if (!VectorContains(so2, so1.at(i))) unique.push_back(so1.at(i));

            // assign the matrix element according to Slater-Condon rules
            if (diff == 0) {
                for (int so : so1) elem += Hms(so, so);
                for (size_t i = 0; i < so1.size(); i++) {
                    for (size_t j = 0; j < so1.size(); j++) {
                        elem += 0.5 * (Jms(so1.at(i), so1.at(i), so1.at(j), so1.at(j)) - Jms(so1.at(i), so1.at(j), so1.at(j), so1.at(i)));
                    }
                }
            } else if (diff == 1) {
                elem += Hms(unique.at(0), unique.at(1));
                for (int so : common) {
                    elem += Jms(unique.at(0), unique.at(1), so, so) - Jms(unique.at(0), so, so, unique.at(1));
                }
            } else if (diff == 2) {
                elem = Jms(unique.at(0), unique.at(2), unique.at(1), unique.at(3)) - Jms(unique.at(0), unique.at(3), unique.at(1), unique.at(2));
            }

            // return element
            Hci(i, j) = std::pow(-1, swaps) * elem;
        }
    }

    Eigen::SelfAdjointEigenSolver<EigenMatrix<>> solver(Hci);

    EigenVector<> Eci = solver.eigenvalues();

    std::cout << solver.eigenvalues().transpose() << std::endl;

    // print the final energy and total time
    std::printf("\nFINAL SINGLE POINT ENERGY: %.14f\n\nTOTAL TIME: %s\n", Eci(0) + system.nuclearRepulsion(), Timer::Format(Timer::Elapsed(start)).c_str());
}
