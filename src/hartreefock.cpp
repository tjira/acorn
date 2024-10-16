#include "hartreefock.h"

torch::Tensor HartreeFock::get_density(const System& system, const torch::Tensor& C_MO) const {
    // define the prefactor and occupied orbitals (spatial or spin)
    double factor = input.generalized ? 1.0 : 2.0; int occupied_orbitals = input.generalized ? system.electrons() : system.electrons() / 2;

    // return the density matrix
    return factor * C_MO.index({"...", Slice(None, occupied_orbitals)}).mm(C_MO.index({"...", Slice(None, occupied_orbitals)}).swapaxes(0, 1));
}

double HartreeFock::get_energy(const torch::Tensor& H_AO, const torch::Tensor& F_AO, const torch::Tensor& D_MO) const {
    return 0.5 * (D_MO * (H_AO + F_AO)).sum().item<double>();
}

torch::Tensor HartreeFock::get_fock(const torch::Tensor& H_AO, const torch::Tensor& J_AO, const torch::Tensor& D_MO) const {
    return H_AO + torch::einsum("ijkl,ij->kl", {J_AO - (input.generalized ? 1.0 : 0.5) * J_AO.swapaxes(1, 3), D_MO});
}

torch::Tensor HartreeFock::gradient(const System& system, const torch::Tensor& dT_AO, const torch::Tensor& dV_AO, const torch::Tensor& dS_AO, const torch::Tensor& dJ_AO, const torch::Tensor& F_AO, const torch::Tensor& C_MO) const {
    // throw an error for generalized calculations
    if (input.generalized) throw std::runtime_error("GRADIENTS NOT SUPPORTED FOR GENERALIZED HARTREE-FOCK");

    // extract the one-integral integral derivative slices
    torch::Tensor dT_AO_1 = dT_AO.index({Slice(), Slice(), Slice(None, 3)}), dV_AO_1 = dV_AO.index({Slice(), Slice(), Slice(None, 3)}), dS_AO_1 = dS_AO.index({Slice(), Slice(), Slice(None, 3)});

    // extract the two-integral integral derivative slices and the density matrix
    torch::Tensor dJ_AO_1 = dJ_AO.index({Slice(), Slice(), Slice(), Slice(), Slice(None, 3)}), D_MO = get_density(system, C_MO);

    // extract the orbital energies, number of occupied orbitals and shell map
    torch::Tensor eps = Transform::SingleSpatial(F_AO, C_MO).diagonal(); int nocc = system.electrons() / 2; auto atom2shell = system.get_shells().atom2shell(system.get_atoms());

    // initialize the gradient matrix and weighted density matrix
    torch::Tensor G = torch::zeros({(int)system.get_atoms().size(), 3}, torch::kDouble), W = torch::zeros_like(C_MO);

    // fill the weighted density matrix
    for (int i = 0; i < W.size(0); i++) for (int j = 0; j < W.size(1); j++) {
        W.index_put_({i, j}, 2 * (C_MO.index({i, Slice(None, nocc)}) * C_MO.index({j, Slice(None, nocc)}) * eps.index({Slice(None, nocc)})).sum());
    }

    // contract the ERI derivative with the second density matrix
    torch::Tensor dERI =  torch::einsum("kl,ijklx->ijx", {D_MO, dJ_AO_1 - 0.5 * dJ_AO_1.swapaxes(1, 3)});

    // calculate the gradient matrix
    for (int i = 0, si = 0, ss = 0; i < G.size(0); i++, si += ss, ss = 0) {

        // calculate number of shells for current atom
        for (long shell : atom2shell.at(i)) ss += system.get_shells().at(shell).size();

        // define the hamiltonian derivative as the nuclear derivative terms
        torch::Tensor dH_AO = dV_AO.index({Slice(), Slice(), Slice(6 + i * 3, 9 + i * 3)});

        // add the orbital derivative terms
        dH_AO.index({Slice(), Slice(si, si + ss), Slice(None, 3)}) += (dT_AO_1 + dV_AO_1).index({Slice(si, si + ss), Slice(), Slice(None, 3)}).swapaxes(0, 1);
        dH_AO.index({Slice(si, si + ss), Slice(), Slice(None, 3)}) += (dT_AO_1 + dV_AO_1).index({Slice(si, si + ss), Slice(), Slice(None, 3)}).swapaxes(0, 0);

        // add the first term of the gradient matrix
        G.index({i, Slice()}) += torch::einsum("ijx,ij->x", {dH_AO, D_MO});

        // add the second and third term of the gradient matrix
        G.index({i, Slice()}) += 2 * torch::einsum("ijx,ij->x", {dERI.index({Slice(si, si + ss), Slice(), Slice()}), D_MO.index({Slice(si, si + ss), Slice()})});
        G.index({i, Slice()}) -= 2 * torch::einsum("ijx,ij->x", {dS_AO_1.index({Slice(si, si + ss), Slice(), Slice()}), W.index({Slice(si, si + ss), Slice()})});
    }

    // return the gradient matrix
    return G;
}


std::tuple<torch::Tensor, torch::Tensor, double> HartreeFock::run(const System& system, const torch::Tensor& H_AO, const torch::Tensor& S_AO, const torch::Tensor& J_AO, torch::Tensor D_MO) const {
    // throw an error if the system is open shell and not specified as generalized
    if (system.get_multi() != 1 && !input.generalized) throw std::runtime_error("OPEN SHELL SYSTEMS REQUIRE GENERALIZED HARTREE-FOCK");

    // define the matrices in the spinorbital basis
    torch::Tensor J_AS, H_AS, S_AS, D_MS;

    // transform the coulomb tensor
    if (input.generalized) J_AS = torch::kron(torch::eye(2, 2, torch::kDouble), torch::kron(torch::eye(2, 2, torch::kDouble), J_AO).swapaxes(0, 3).swapaxes(1, 2).contiguous());

    // transform the one-electron integrals
    if (input.generalized) H_AS = torch::kron(torch::eye(2, 2, torch::kDouble), H_AO);
    if (input.generalized) S_AS = torch::kron(torch::eye(2, 2, torch::kDouble), S_AO);

    // transform the density matrix
    if (input.generalized) D_MS = torch::kron(torch::eye(2, 2, torch::kDouble), D_MO);

    // run the Hartree-Fock calculation
    return input.generalized ? scf(system, H_AS, S_AS, J_AS, D_MS) : scf(system, H_AO, S_AO, J_AO, D_MO);
}

std::tuple<torch::Tensor, torch::Tensor, double> HartreeFock::scf(const System& system, const torch::Tensor& H, const torch::Tensor& S, const torch::Tensor& J, torch::Tensor D) const {
    // define the Fock matrix, orbital coefficients, electron energy and containers for errors and Fock matrices
    torch::Tensor F = H, C = torch::zeros_like(D); double energy_hf = 0; std::vector<torch::Tensor> errors, focks;

    // calculate the complementary X matrix as square root of the inverse of the overlap matrix
    auto [SVAL, SVEC] = torch::linalg::eigh(S, "L"); torch::Tensor X = SVEC.mm(torch::diag(1 / torch::sqrt(SVAL))).mm(SVEC.swapaxes(0, 1));

    // print the header
    std::printf("%s HF ENERGY CALCULATION\n%6s %20s %8s %8s %12s\n", input.generalized ? "GENERALIZED" : "RESTRICTED", "ITER", "ENERGY", "|dE|", "|dD|", "TIMER");

    // start the SCF procedure
    for (int i = 0; i < input.max_iter; i++) {

        // define the iteration timer
        Timepoint iteration_timer = Timer::Now();

        // calculate the Fock matrix and save the previous values
        F = get_fock(H, J, D); torch::Tensor D_P = D; double energy_hf_prev = energy_hf;

        // calculate the error vector and append it to the container along with the Fock matrix
        torch::Tensor error = S.mm(D).mm(F) - F.mm(D).mm(S); if (i) errors.push_back(error), focks.push_back(F);

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
            F = c.index({0}) * focks.at(0); for (int j = 1; j < input.diis_size; j++) F += c.index({j}) * focks.at(j);
        }

        // solve the Roothaan equations
        torch::Tensor E; std::tie(E, C) = torch::linalg::eigh(X.mm(F).mm(X), "L"); C = X.mm(C);

        // calculate the new density and energy
        D = get_density(system, C), energy_hf = get_energy(H, F, D);

        // calculate the errors
        double error_energy = std::abs(energy_hf - energy_hf_prev), error_density = (D - D_P).norm().item<double>();

        // print the iteration info
        std::printf("%6d %20.14f %.2e %.2e %s %s\n", i + 1, energy_hf, error_energy, error_density, Timer::Format(Timer::Elapsed(iteration_timer)).c_str(), input.diis_size && i >= input.diis_size ? "DIIS" : "");

        // finish if covergence reached
        if (std::abs(energy_hf - energy_hf_prev) < input.threshold && (D - D_P).norm().item<double>() < input.threshold) break;

        // throw an exception if the maximum number of iterations reached
        if (i == input.max_iter - 1) throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN THE SCF REACHED");
    }

    // return the converged results
    return {F, C, energy_hf};
}
