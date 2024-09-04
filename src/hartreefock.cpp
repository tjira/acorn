#include "hartreefock.h"

torch::Tensor HartreeFock::get_density(const System& system, const torch::Tensor& C_MO) {
    return 2 * C_MO.index({"...", Slice(None, system.occupied_spatial_orbitals())}).mm(C_MO.index({"...", Slice(None, system.occupied_spatial_orbitals())}).swapaxes(0, 1));
}

double HartreeFock::get_energy(const torch::Tensor& H_AO, const torch::Tensor& F_AO, const torch::Tensor& D_MO) {
    return 0.5 * (D_MO * (H_AO + F_AO)).sum().item<double>();
}

torch::Tensor HartreeFock::get_fock(const torch::Tensor& H_AO, const torch::Tensor& J_AO, const torch::Tensor& D_MO) {
    return H_AO + torch::einsum("ijkl,ij->kl", {J_AO - 0.5 * J_AO.swapaxes(0, 3), D_MO});
}

torch::Tensor HartreeFock::run(const System& system, const torch::Tensor& H_AO, const torch::Tensor& S_AO, const torch::Tensor& J_AO, torch::Tensor D_MO) const {
    // define the Fock matrix, orbital coefficients, electron energy and containers for errors and Fock matrices
    torch::Tensor F_AO = H_AO, C_MO = torch::zeros_like(D_MO); double energy_hf = 0; std::vector<torch::Tensor> errors, focks;

    // calculate the complementary X matrix as square root of the inverse of the overlap matrix
    auto [SVAL, SVEC] = torch::linalg::eigh(S_AO, "L"); torch::Tensor X = SVEC.mm(torch::diag(1 / torch::sqrt(SVAL))).mm(SVEC.swapaxes(0, 1));

    // print the header
    std::printf("HF ENERGY CALCULATION\n%6s %20s %8s %8s %12s\n", "ITER", "ENERGY", "|dE|", "|dD|", "TIMER");

    // start the SCF procedure
    for (int i = 0; i < input.max_iter; i++) {

        // define the iteration timer
        Timepoint iteration_timer = Timer::Now();

        // calculate the Fock matrix and save the previous values
        F_AO = get_fock(H_AO, J_AO, D_MO); torch::Tensor D_MO_P = D_MO; double energy_hf_prev = energy_hf;

        // calculate the error vector and append it to the container along with the Fock matrix
        torch::Tensor error = S_AO.mm(D_MO).mm(F_AO) - F_AO.mm(D_MO).mm(S_AO); if (i) errors.push_back(error), focks.push_back(F_AO);

        // truncate the error and Fock vector containers if they exceed the DIIS size
        if (i > input.diis_size) errors.erase(errors.begin()), focks.erase(focks.begin());

        // perform DIIS extrapolation
        if (input.diis_size && i >= input.diis_size) {

            // define the DIIS subspace matrices
            torch::Tensor B = torch::ones ({input.diis_size + 1, input.diis_size + 1}, torch::kDouble);
            torch::Tensor b = torch::zeros({input.diis_size + 1                     }, torch::kDouble);

            // fill the edge values DIIS subspace matrices
            B.index({input.diis_size, input.diis_size}) = 0, b.index({input.diis_size}) = 1;

            // fill the DIIS matrix
            for (int j = 0; j < input.diis_size; j++) for (int k = j; k < input.diis_size; k++) B.index({j, k}) = (errors.at(j) * errors.at(k)).sum().item<double>(), B.index({k, j}) = B.index({j, k});

            // solve the DIIS equations for the coefficients
            torch::Tensor c = torch::linalg::solve(B, b, true);

            // extrapolate the Fock matrix
            F_AO = c.index({0}) * focks.at(0); for (int j = 1; j < input.diis_size; j++) F_AO += c.index({j}) * focks.at(j);
        }

        // solve the Roothaan equations
        torch::Tensor E_MO; std::tie(E_MO, C_MO) = torch::linalg::eigh(X.mm(F_AO).mm(X), "L"); C_MO = X.mm(C_MO);

        // calculate the new density and energy
        D_MO = get_density(system, C_MO), energy_hf = get_energy(H_AO, F_AO, D_MO);

        // calculate the errors
        double error_energy = std::abs(energy_hf - energy_hf_prev), error_density = (D_MO - D_MO_P).norm().item<double>();

        // print the iteration info
        std::printf("%6d %20.14f %.2e %.2e %s %s\n", i + 1, energy_hf, error_energy, error_density, Timer::Format(Timer::Elapsed(iteration_timer)).c_str(), input.diis_size && i >= input.diis_size ? "DIIS" : "");

        // finish if covergence reached
        if (std::abs(energy_hf - energy_hf_prev) < input.threshold && (D_MO - D_MO_P).norm().item<double>() < input.threshold) break;

        // throw an exception if the maximum number of iterations reached
        if (i == input.max_iter - 1) throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN THE SCF REACHED");
    }

    // return the converged orbital coefficients
    return C_MO;
}