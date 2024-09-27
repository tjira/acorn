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

std::vector<std::vector<int>> ConfigurationInteraction::generate_all_generalized_configurations(const System& system) const {
    // define the number of electrons and orbitals
    int electrons = input.cas.empty() ? system.electrons() : input.cas.at(0), orbitals = input.cas.empty() ? system.basis_functions() : input.cas.at(1);

    // check for nonsense active space
    if (electrons > system.electrons() || orbitals > system.basis_functions() || orbitals > electrons / 2 + system.virtual_spinorbitals() / 2) throw std::invalid_argument("YOUR ACTIVE SPACE IS NONSENSE");

    // define the fixed spinorbitals
    std::vector<int> fixed(system.electrons() - electrons); std::iota(fixed.begin(), fixed.end(), 0);

    // generate the configurations in the active space
    std::vector<std::vector<int>> combinations = Combinations(2 * orbitals, electrons);

    // edit the configurations
    for (size_t i = 0; i < combinations.size(); i++) {

        // since the Combinations function generates from 0, we need to add the number of fixed electrons to each element
        std::transform(combinations.at(i).begin(), combinations.at(i).end(), combinations.at(i).begin(), [&](int x) {return x + system.electrons() - electrons;});

        // insert the fixed spinorbitals at the beginning of the vector
        combinations.at(i).insert(combinations.at(i).begin(), fixed.begin(), fixed.end());
    }

    return combinations;
}

std::vector<std::vector<int>> ConfigurationInteraction::generate_all_restricted_configurations(const System& system) const {
    // define the number of electrons and orbitals
    int electrons = input.cas.empty() ? system.electrons() : input.cas.at(0), orbitals = input.cas.empty() ? system.basis_functions() : input.cas.at(1);

    // check for nonsense active space
    if (electrons > system.electrons() || orbitals > system.basis_functions() || orbitals > electrons / 2 + system.virtual_spinorbitals() / 2) throw std::invalid_argument("YOUR ACTIVE SPACE IS NONSENSE");

    // define the determinant vector and single spin electron configurations
    std::vector<std::vector<int>> configs = Combinations(orbitals, electrons / 2); std::vector<std::vector<int>> combinations(configs.size() * configs.size());

    // define the fixed spinorbitals and the orbital offset
    std::vector<int> fixed(system.electrons() - electrons); std::iota(fixed.begin(), fixed.end(), 0);

    // generate all possible combinations of alpha and beta electrons
    for (size_t i = 0; i < configs.size(); i++) for (size_t j = 0; j < configs.size(); j++) {
        for (int k = 0; k < electrons / 2; k++) combinations.at(i * configs.size() + j).push_back(2 * configs.at(i).at(k) + 0 + system.electrons() - electrons);
        for (int k = 0; k < electrons / 2; k++) combinations.at(i * configs.size() + j).push_back(2 * configs.at(j).at(k) + 1 + system.electrons() - electrons);
    }

    // sort and add the fixed spinorbitals
    for (size_t i = 0; i < combinations.size(); i++) std::sort(combinations.at(i).begin(), combinations.at(i).end()), combinations.at(i).insert(combinations.at(i).begin(), fixed.begin(), fixed.end());

    // return the determinant vector
    return combinations;
}

std::string ConfigurationInteraction::get_name() const {
    // initialize the basi name
    std::string name = "CI";

    // add the F if FCI performed
    if (input.cas.size() == 0) name = "F" + name;

    // add the number of electrons and orbitals
    if (input.cas.size() == 2) name = "CAS" + name + "(" + std::to_string(input.cas.at(0)) + ", " + std::to_string(input.cas.at(1)) + ")";

    // return the name
    return name;
}

std::tuple<torch::Tensor, torch::Tensor> ConfigurationInteraction::run(const System& system, const torch::Tensor& H_MS, const torch::Tensor& J_MS_AP) const {
    // start the timer for determinant generation
    Timepoint generation_timer = Timer::Now(); std::printf("\nGENERATING ALL CONFIGURATIONS: "); std::flush(std::cout);

    // generate the configurations
    std::vector<std::vector<int>> dets = input.triplet ? generate_all_generalized_configurations(system) : generate_all_restricted_configurations(system);

    // initialize the data for the hamiltonian
    std::vector<double> hamiltonian_data(dets.size() * dets.size(), 0);

    // print the time taken to generate the configurations
    std::printf("%s\n\nNUMBER OF DETERMINANTS GENERATED: %lu\n\nGENERATING THE CI HAMILTONIAN: ", Timer::Format(Timer::Elapsed(generation_timer)).c_str(), dets.size()); std::flush(std::cout);

    // create element accesors for the Hamiltonian and J tensor
    torch::TensorAccessor<double, 2> H_MS_accessor = H_MS.accessor<double, 2>(); torch::TensorAccessor<double, 4> J_MS_AP_accessor = J_MS_AP.accessor<double, 4>();

    // start the timer for Hamiltonian filling
    Timepoint hamiltonian_timer = Timer::Now();

    // fill the CI Hamiltonian
    #pragma omp parallel for num_threads(nthread)
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
    std::printf("%s\nDIAGONALIZING THE HAMILTONIAN: ", Timer::Format(Timer::Elapsed(hamiltonian_timer)).c_str()); std::flush(std::cout);

    // start the timer for diagonalization
    Timepoint diagonalization_timer = Timer::Now();

    // diagonalize the CI Hamiltonian
    auto [E_CI, C_CI] = torch::linalg::eigh(torch::from_blob(hamiltonian_data.data(), {(int)dets.size(), (int)dets.size()}, torch::kDouble), "L");

    // print the time taken to diagonalize the Hamiltonian
    std::printf("%s\n", Timer::Format(Timer::Elapsed(diagonalization_timer)).c_str());

    // return the results
    return {E_CI, C_CI};
}
