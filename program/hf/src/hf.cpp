#include "hf.h"

void Acorn::HF::run(const Options& opt, std::vector<timepoint>& timers) {
    // start the timer for integral loading
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of integral loading
    std::cout << "SYSTEM AND INTEGRALS IN AO BASIS READING: " << std::flush;

    // load the integrals in AO basis
    torch::Tensor V = torch::ReadTensor("V_AO.mat");
    torch::Tensor T = torch::ReadTensor("T_AO.mat");
    torch::Tensor S = torch::ReadTensor("S_AO.mat");
    torch::Tensor J = torch::ReadTensor("J_AO.mat");

    // print the time for integral loading
    std::cout << eltime(timers.at(1)) << std::endl;

    // Hamiltonian, Fock matrix energy and DIIS containers
    torch::Tensor H = T + V, F = H; double Eel = 0; std::vector<torch::Tensor> es, fs;

    // initial guess for the density and MO coefficients
    torch::Tensor Cmo = torch::zeros({S.sizes().at(0), S.sizes().at(1)}, torch::kDouble);
    torch::Tensor Dmo = torch::zeros({S.sizes().at(0), S.sizes().at(1)}, torch::kDouble);
    torch::Tensor Emo = torch::zeros({S.sizes().at(0), 1              }, torch::kDouble);

    // calculate the complementary X matrix as square root of the inverse of the overlap matrix
    auto[SVAL, SVEC] = torch::linalg::eigh(S, "L"); torch::Tensor X = SVEC.mm(torch::diag(1 / torch::sqrt(SVAL))).mm(SVEC.swapaxes(0, 1));

    // print the header
    std::printf("\n%6s %20s %8s %8s %12s\n", "ITER", "ENERGY", "|dE|", "|dD|", "TIME");

    // start the SCF procedure
    for (int i = 0; i < opt.iters; i++) {

        // start the timer
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // calculate the Fock matrix and define previous values
        F = H + torch::einsum("ijkl,ij->kl", {J - 0.5 * J.swapaxes(0, 3), Dmo}); torch::Tensor Dmop = Dmo; double Eelp = Eel;

        // calculate the error vector and append it to the container along with the Fock matrix
        torch::Tensor e = S.mm(Dmo).mm(F) - F.mm(Dmo).mm(S); if (i) es.push_back(e), fs.push_back(F);

        // truncate the error and Fock vector containers
        if (i > opt.diis) {es.erase(es.begin()), fs.erase(fs.begin());}

        // perform DIIS extrapolation
        if (opt.diis && i >= opt.diis) {

            // define the DIIS subspace matrices
            torch::Tensor B = torch::ones( {opt.diis + 1, opt.diis + 1}, torch::kDouble);
            torch::Tensor b = torch::zeros({opt.diis + 1              }, torch::kDouble);

            // fill the edge values DIIS subspace matrices
            B.index({opt.diis, opt.diis}) = 0, b.index({opt.diis}) = 1;

            // fill the DIIS matrix
            for (int j = 0; j < opt.diis; j++) for (int k = j; k < opt.diis; k++) B.index({j, k}) = (es.at(j) * es.at(k)).sum().item<double>(), B.index({k, j}) = B.index({j, k});

            // solve the DIIS equations
            torch::Tensor c = torch::linalg::solve(B, b, true);

            // extrapolate the Fock matrix
            F = c.index({0}) * fs.at(0); for (int j = 1; j < opt.diis; j++) F += c.index({j}) * fs.at(j);
        }

        // solve the Roothaan equations
        std::tie(Emo, Cmo) = torch::linalg::eigh(X.mm(F).mm(X), "L"); Cmo = X.mm(Cmo);

        // calculate the new density
        Dmo = 2 * Cmo.index({"...", Slice(None, opt.system.nocc())}).mm(Cmo.index({"...", Slice(None, opt.system.nocc())}).swapaxes(0, 1));

        // calculate the new energy
        Eel = 0.5 * (Dmo * (H + F)).sum().item<double>();

        // print the iteration info
        std::printf("%6d %20.14f %.2e %.2e %s %s\n", i + 1, Eel, std::abs(Eel - Eelp), (Dmo - Dmop).norm().item<double>(), eltime(timers.at(1)).c_str(), opt.diis && i >= opt.diis ? "DIIS" : "");

        // finish if covergence reached
        if (std::abs(Eel - Eelp) < opt.thresh && (Dmo - Dmop).norm().item<double>() < opt.thresh) {std::cout << std::endl; break;}
        else if (i == opt.iters - 1) throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN THE SCF REACHED");
    }

    // calculate the nuclear repulsion
    torch::Tensor N = torch::zeros({2}); N.index({0}) = opt.system.nocc(), N.index({1}) = opt.system.nuclearRepulsion();

    // start the timer for writing the results
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of the writing of final results
    std::cout << "HARTREE-FOCK RESULTS WRITING: " << std::flush;

    // save the final matrices
    torch::WriteTensor("C_MO.mat", Cmo);
    torch::WriteTensor("E_MO.mat", Emo);
    torch::WriteTensor("D_MO.mat", Dmo);
    torch::WriteTensor("F_AO.mat", F  );
    torch::WriteTensor("N.mat",    N  );

    // print the time for writing the results
    std::cout << eltime(timers.at(1)) << std::endl;

    // print the final energy and total time
    std::printf("\nFINAL SINGLE POINT ENERGY: %.14f\n\nTOTAL TIME: %s\n", Eel + N.index({1}).item<double>(), eltime(timers.at(0)).c_str());
}
