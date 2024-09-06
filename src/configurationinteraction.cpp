#include "configurationinteraction.h"

std::tuple<std::vector<int>, std::vector<int>, int> ConfigurationInteraction::align_determinants(std::vector<int> deta, const std::vector<int>& detb) {
    // define the number of swaps
    int swaps = 0;

    // align the determinants
    for (size_t i = 0; i < deta.size(); i++) if (deta.at(i) != detb.at(i)) {
        for (size_t j = 0; j < deta.size(); j++) if (deta.at(i) == detb.at(j)) {
            std::swap(deta.at(i), deta.at(j)), swaps++;
        }
    }

    // return the aligned determinants and the number of swaps
    return {deta, detb, swaps};
}

std::vector<std::vector<int>> ConfigurationInteraction::generate_all_configurations(const System& system) {
    // define the determinant vector and single electron configurations
    std::vector<std::vector<int>> configs = Combinations(system.basis_functions(), system.occupied_spatial_orbitals()); std::vector<std::vector<int>> dets(configs.size() * configs.size());

    // generate all possible combinations of alpha and beta electrons
    for (size_t i = 0; i < configs.size(); i++) for (size_t j = 0; j < configs.size(); j++) {
        for (int k = 0; k < system.occupied_spatial_orbitals(); k++) dets.at(i * configs.size() + j).push_back(2 * configs.at(i).at(k) + 0);
        for (int k = 0; k < system.occupied_spatial_orbitals(); k++) dets.at(i * configs.size() + j).push_back(2 * configs.at(j).at(k) + 1);
    }

    // return the determinant vector
    return dets;
}

std::tuple<std::vector<int>, std::vector<int>> ConfigurationInteraction::get_common_and_unique_spinorbitals(std::vector<int> deta, const std::vector<int>& detb) {
    // define the determinants
    std::vector<int> common, unique;

    // fill the common spinorbitals and calculate the number of differences
    for (size_t i = 0; i < deta.size(); i++) if (deta.at(i) == detb.at(i)) common.push_back(deta.at(i));

    // fill the unique spinorbitals vector
    for (size_t i = 0; i < deta.size(); i++) if (std::find(deta.begin(), deta.end(), detb.at(i)) == deta.end()) unique.push_back(detb.at(i));
    for (size_t i = 0; i < detb.size(); i++) if (std::find(detb.begin(), detb.end(), deta.at(i)) == detb.end()) unique.push_back(deta.at(i));

    // return the common and unique spinorbitals
    return {common, unique};
}

double ConfigurationInteraction::slater_condon_rules(std::vector<int> deta, const std::vector<int>& detb, const torch::TensorAccessor<double, 2>& H_MS_accessor, const torch::TensorAccessor<double, 4>& J_MS_AP_accessor) {
    // define the element
    double elem = 0;

    // get the common and unique spinorbitals
    auto [common, unique] = get_common_and_unique_spinorbitals(deta, detb);

    // apply the Slater-Condon rules to the matrix elemetnt
    if (unique.size() == 0) {
        for (int so : deta) elem += H_MS_accessor[so][so];
        for (size_t i = 0; i < deta.size(); i++) for (size_t j = 0; j < deta.size(); j++) {
            elem += 0.5 * J_MS_AP_accessor[deta.at(i)][deta.at(j)][deta.at(i)][deta.at(j)];
        }
    } else if (unique.size() == 2) {
        elem += H_MS_accessor[unique.at(0)][unique.at(1)];
        for (int so : common) elem += J_MS_AP_accessor[unique.at(0)][so][unique.at(1)][so];
    } else if (unique.size() == 4) {
        elem = J_MS_AP_accessor[unique.at(0)][unique.at(1)][unique.at(2)][unique.at(3)];
    }

    // return the matrix element
    return elem;
}

std::string ConfigurationInteraction::get_name() const {
    // initialize the basi name
    std::string name = "CI";

    // add the excitation level to the name
    if (input.excitation.size() == 0) name = "F" + name;

    // return the name
    return name;
}

std::tuple<torch::Tensor, torch::Tensor> ConfigurationInteraction::run(const System& system, const torch::Tensor& H_MS, const torch::Tensor& J_MS_AP) const {
    // throw an error if the excitation level is not supported
    if (input.excitation.size() != 0) throw std::runtime_error("PROVIDED CI EXCITATION LEVEL NOT SUPPORTED");

    // start the timer for determinant generation
    Timepoint generation_timer = Timer::Now(); std::printf("\nGENERATING ALL CONFIGURATIONS: ");

    // generate the configurations and define the CI hamiltonian data
    std::vector<std::vector<int>> dets = generate_all_configurations(system); std::vector<double> hamiltonian_data(dets.size() * dets.size(), 0);

    // print the time taken to generate the configurations
    std::printf("%s\nGENERATING THE CI HAMILTONIAN: ", Timer::Format(Timer::Elapsed(generation_timer)).c_str());

    // create element accesors for the Hamiltonian and J tensor
    torch::TensorAccessor<double, 2> H_MS_accessor = H_MS.accessor<double, 2>(); torch::TensorAccessor<double, 4> J_MS_AP_accessor = J_MS_AP.accessor<double, 4>();

    // start the timer for Hamiltonian filling
    Timepoint hamiltonian_timer = Timer::Now();

    // fill the CI Hamiltonian
    for (size_t i = 0; i < dets.size(); i++) {
        for (size_t j = i; j < dets.size(); j++) {
            
            // align the determinants
            auto [deta, detb, swaps] = align_determinants(dets.at(i), dets.at(j));

            // assign the element multiplied by the correct sign
            hamiltonian_data.at(i * dets.size() + j) = std::pow(-1, swaps) * slater_condon_rules(deta, detb, H_MS_accessor, J_MS_AP_accessor);

            // symmetrize the matrix
            hamiltonian_data.at(j * dets.size() + i) = hamiltonian_data.at(i * dets.size() + j);
        }
    }

    // print the time taken to fill the Hamiltonian
    std::printf("%s\nDIAGONALIZING THE HAMILTONIAN: ", Timer::Format(Timer::Elapsed(hamiltonian_timer)).c_str());

    // start the timer for diagonalization
    Timepoint diagonalization_timer = Timer::Now();

    // diagonalize the CI Hamiltonian
    auto [E_CI, C_CI] = torch::linalg::eigh(torch::from_blob(hamiltonian_data.data(), {(int)dets.size(), (int)dets.size()}, torch::kDouble), "L");

    // print the time taken to diagonalize the Hamiltonian
    std::printf("%s\n", Timer::Format(Timer::Elapsed(diagonalization_timer)).c_str());

    // return the results
    return {E_CI, C_CI};
}
