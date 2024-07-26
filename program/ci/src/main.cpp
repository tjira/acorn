#include "linalg.h"
#include "timer.h"
#include <argparse.hpp>

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

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // load the integrals in AO basis and system from disk
    MEASURE("SYSTEM AND INTEGRALS IN MS BASIS READING: ",
        Matrix    Vms = Eigen::LoadMatrix("V_MS.mat");
        Matrix    Tms = Eigen::LoadMatrix("T_MS.mat");
        Tensor<4> Jms = Eigen::LoadTensor("J_MS.mat");
        Vector    N   = Eigen::LoadMatrix("N.mat"   );
    )

    // define the number of basis functions and occupied orbitals, Hamiltonian in MS basis and determinant container
    int nbf = Vms.rows() / 2, nocc = N(0); Matrix Hms = Tms + Vms; std::vector<std::vector<int>> dets;

    MEASURE("GENERATING ALL POSSIBLE DETERMINANTS:     ",
        for (const std::vector<int>& alpha : Combinations(nbf, nocc)) {
            for (const std::vector<int>& beta : Combinations(nbf, nocc)) {
                dets.push_back({});

                for (int i = 0; i < nocc; i++) dets.back().push_back(2 * alpha.at(i));
                for (int i = 0; i < nocc; i++) dets.back().push_back(2 * beta.at(i) + 1);
            }
        }
    ) std::cout << "\nNUMBER OF DETERMINANTS GENERATED: " << dets.size() << std::endl;

    // print the label of Hamiltonian filling timer
    std::cout << "\nFILLING THE CI HAMILTONIAN:  " << std::flush; tp = Timer::Now();

    // define the CI Hamiltonian
    Matrix Hci = Matrix::Zero(dets.size(), dets.size());

    // fill the CI Hamiltonian
    for (int i = 0; i < Hci.rows(); i++) {
        for (int j = i; j < Hci.cols(); j++) {
            
            // extract the determinants, define the element, differences and swaps and initialize the common and unique spinorbitals
            std::vector<int> deta = dets.at(i), detb = dets.at(j); int diff = 0, swaps = 0; double elem = 0; std::vector<int> common, unique;

            // align the determinants
            for (size_t i = 0; i < deta.size(); i++) {
                if (deta.at(i) != detb.at(i)) {
                    for (size_t j = 0; j < deta.size(); j++) {
                        if (deta.at(i) == detb.at(j)) {
                            std::swap(deta.at(i), deta.at(j)), swaps++;
                        }
                    }
                }
            }

            // fill the common spinorbitals and calculate the number of differences
            for (int i = 0; i < 2 * nocc; i++) {if (deta.at(i) == detb.at(i)) common.push_back(deta.at(i)); else diff++;}

            // fill the unique spinorbitals vector
            for (int i = 0; i < 2 * nocc; i++) if (std::find(deta.begin(), deta.end(), detb.at(i)) == deta.end()) unique.push_back(detb.at(i));
            for (int i = 0; i < 2 * nocc; i++) if (std::find(detb.begin(), detb.end(), deta.at(i)) == detb.end()) unique.push_back(deta.at(i));

            // apply the Slater-Condon rules to the matrix elemetnt
            if (diff == 0) {
                for (int so : deta) elem += Hms(so, so);
                for (size_t i = 0; i < deta.size(); i++) {
                    for (size_t j = 0; j < deta.size(); j++) {
                        elem += 0.5 * (Jms(deta.at(i), deta.at(i), deta.at(j), deta.at(j)) - Jms(deta.at(i), deta.at(j), deta.at(j), deta.at(i)));
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

            // return element multiplied by the correct sign
            Hci(i, j) = std::pow(-1, swaps) * elem; Hci(j, i) = Hci(i, j);
        }
    }

    // print the time taken to fill the CI Hamiltonian
    std::cout << Timer::Format(Timer::Elapsed(tp)) << std::endl;

    // find the eigenvalues and eigenvectors of the CI Hamiltonian
    MEASURE("SOLVING THE CI EIGENPROBLEM: ", Eigen::SelfAdjointEigenSolver<Matrix> solver(Hci)) Matrix Cci = solver.eigenvectors(), Eci = solver.eigenvalues();

    // write the results
    MEASURE("WRITING THE RESULT MATRICES: ",
        Eigen::Write("H_CI.mat", Hci);
        Eigen::Write("C_CI.mat", Cci);
        Eigen::Write("E_CI.mat", Eci);
    )

    // print the final energy and total time
    std::printf("\nFINAL SINGLE POINT ENERGY: %.14f\n\nTOTAL TIME: %s\n", Eci(0) + N(1), Timer::Format(Timer::Elapsed(start)).c_str());
}
