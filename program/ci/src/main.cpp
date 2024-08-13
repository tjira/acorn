#include "tensor.h"
#include <argparse.hpp>

#define FORMAT(T) [&](long ms) {char s[99]; std::sprintf(s, "%02ld:%02ld:%02ld.%03ld", ms / 3600000, ms % 3600000 / 60000, ms % 60000 / 1000, ms % 1000); return std::string(s);}(T)

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
    argparse::ArgumentParser program("Acorn Configuration Interaction Program", "1.0", argparse::default_arguments::none);

    // define the timers
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> timers(2, std::chrono::high_resolution_clock().now()); auto tp = timers.at(0);
    auto elapsed = [](auto ms) {return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock().now() - ms).count();};

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // start the timer for the integral loading
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of integral loading
    std::cout << "SYSTEM AND INTEGRALS IN MS BASIS READING: " << std::flush;

    // load the integrals in AO basis and system from disk
    torch::Tensor Vms = torch::ReadTensor("V_MS.mat");
    torch::Tensor Tms = torch::ReadTensor("T_MS.mat");
    torch::Tensor Jms = torch::ReadTensor("J_MS.mat");
    torch::Tensor N   = torch::ReadTensor("N.mat"   );

    // print the time for integral loading
    std::cout << FORMAT(elapsed(timers.at(1))) << std::endl;

    // define the number of basis functions and occupied orbitals, Hamiltonian in MS basis and determinant container
    int nbf = Vms.sizes().at(0) / 2, nocc = N.index({0}).item<double>(); torch::Tensor Hms = Tms + Vms; std::vector<std::vector<int>> dets;

    // start the timer for generating all possible determinants
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of determinant generation
    std::cout << "GENERATING ALL POSSIBLE DETERMINANTS:     " << std::flush;

    for (const std::vector<int>& alpha : Combinations(nbf, nocc)) {
        for (const std::vector<int>& beta : Combinations(nbf, nocc)) {
            dets.push_back({}); for (int i = 0; i < nocc; i++) {dets.back().push_back(2 * alpha.at(i));} for (int i = 0; i < nocc; i++) {dets.back().push_back(2 * beta.at(i) + 1);}
        }
    } 

    // print the time taken to generate all possible determinants
    std::cout << FORMAT(elapsed(timers.at(1))) << std::endl;

    // print the number of determinants generated
    std::cout << "\nNUMBER OF DETERMINANTS GENERATED: " << dets.size() << std::endl;

    // print the label of Hamiltonian filling timer
    std::cout << "\nFILLING THE CI HAMILTONIAN:  " << std::flush; timers.at(1) = std::chrono::high_resolution_clock().now();

    // define the CI Hamiltonian
    std::vector<double> Hcid(dets.size() * dets.size(), 0);

    // fill the CI Hamiltonian
    for (size_t i = 0; i < dets.size(); i++) {
        for (size_t j = i; j < dets.size(); j++) {
            
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
                for (int so : deta) elem += Hms.index({so, so}).item<double>();
                for (size_t i = 0; i < deta.size(); i++) {
                    for (size_t j = 0; j < deta.size(); j++) {
                        elem += 0.5 * (Jms.index({deta.at(i), deta.at(i), deta.at(j), deta.at(j)}).item<double>() - Jms.index({deta.at(i), deta.at(j), deta.at(j), deta.at(i)})).item<double>();
                    }
                }
            } else if (diff == 1) {
                elem += Hms.index({unique.at(0), unique.at(1)}).item<double>();
                for (int so : common) {
                    elem += Jms.index({unique.at(0), unique.at(1), so, so}).item<double>() - Jms.index({unique.at(0), so, so, unique.at(1)}).item<double>();
                }
            } else if (diff == 2) {
                elem = Jms.index({unique.at(0), unique.at(2), unique.at(1), unique.at(3)}).item<double>() - Jms.index({unique.at(0), unique.at(3), unique.at(1), unique.at(2)}).item<double>();
            }

            // return element multiplied by the correct sign
            Hcid.at(i * dets.size() + j) = std::pow(-1, swaps) * elem; Hcid.at(j * dets.size() + i) = Hcid.at(i * dets.size() + j);
        }
    }

    // convert the CI Hamiltonian to the tensor form
    torch::Tensor Hci = torch::from_blob(Hcid.data(), {(int)dets.size(), (int)dets.size()}, torch::kDouble);

    // print the time taken to fill the CI Hamiltonian
    std::cout << FORMAT(elapsed(timers.at(1))) << std::endl;

    // start the timer for solving the CI eigenproblem
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the label of the CI eigenproblem solver
    std::cout << "SOLVING THE CI EIGENPROBLEM: " << std::flush;

    // find the eigenvalues and eigenvectors of the CI Hamiltonian
    auto[Eci, Cci] = torch::linalg::eigh(Hci, "L");

    // print the time taken to solve the CI eigenproblem
    std::cout << FORMAT(elapsed(timers.at(1))) << std::endl;

    // start the timer for writing the results
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of the writing of final results
    std::cout << "WRITING THE RESULT MATRICES: " << std::flush;

    // save the final matrices
    torch::WriteTensor("H_CI.mat", Hci);
    torch::WriteTensor("C_CI.mat", Cci);
    torch::WriteTensor("E_CI.mat", Eci);

    // print the time for writing the results
    std::cout << FORMAT(elapsed(timers.at(1))) << std::endl;

    // print the final energy and total time
    std::printf("\nFINAL SINGLE POINT ENERGY: %.14f\n\nTOTAL TIME: %s\n", Eci.index({0}).item<double>() + N.index({1}).item<double>(), FORMAT(elapsed(timers.at(0))).c_str());
}
