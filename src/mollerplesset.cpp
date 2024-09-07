#include "mollerplesset.h"

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

double MollerPlesset::evaluate_contraction(const System& system, const std::string& contraction_string, const torch::Tensor& F_MS, const torch::Tensor& J_MS_AP, int order) {
    // define the contribution and split it
    std::vector<std::string> contractions = Split(contraction_string, ';');

    // replace the plus sign with a comma
    std::replace(contractions.at(0).begin(), contractions.at(0).end(), '+', ',');
    std::replace(contractions.at(1).begin(), contractions.at(1).end(), '+', ',');

    // define the tensor slices for the current contributions with axes, array of orbital energy slices and vector of ones for reshaping
    std::vector<torch::Tensor> views; std::vector<torch::Tensor> epsts; std::vector<int64_t> ones = {-1};

    // extract the contraction indices from the unoptimized contraction
    std::vector<std::string> unoptimized_contraction = Split(contractions.at(0), ',');

    // loop over all coulomb contractions
    for (int k = 0; k < order; k++) {

        // define the vector of current coulomb slices
        std::vector<Slice> axes(4); 

        // loop over the axes and assign the slices
        for (int l = 0; l < 4; l++) {
            if (std::string(OCCUPIED_ORBITALS).find(unoptimized_contraction.at(k).at(l)) != std::string::npos) axes.at(l) = Slice(None, system.electrons());
            if (std::string(VIRTUAL_ORBITALS ).find(unoptimized_contraction.at(k).at(l)) != std::string::npos) axes.at(l) = Slice(system.electrons(), None);
        }

        // add the coulomb integral slice to the views
        views.push_back(J_MS_AP.index({axes.at(0), axes.at(1), axes.at(2), axes.at(3)}));
    }

    // loop over energy contractions and add the energy denominators
    for (size_t k = order; k < unoptimized_contraction.size(); k++) {
        views.push_back(Transform::ExcitationEnergyFraction(F_MS, Slice(None, system.electrons()), Slice(system.electrons(), None), unoptimized_contraction.at(k).size() / 2));
    }

    // define the contraction path
    std::vector<std::tuple<std::string, std::array<int, 2>>> path;

    // split the nodes and add them to the path
    for (const std::string& node: Split(contractions.at(1), ':')) {
        path.push_back({Split(node, '/').at(1), {std::stoi(Split(Split(node, '/').at(0), '-').at(0)), std::stoi(Split(Split(node, '/').at(0), '-').at(1))}});
    }

    // contract over every node
    for (const auto& node : path) {

        // contract the nodes
        views.push_back(torch::einsum(std::get<0>(node), {views.at(std::get<1>(node).at(0)), views.at(std::get<1>(node).at(1))}));

        // erase the contracted nodes
        views.erase(views.begin() + std::get<1>(node).at(std::get<1>(node).at(0) < std::get<1>(node).at(1)));
        views.erase(views.begin() + std::get<1>(node).at(std::get<1>(node).at(0) > std::get<1>(node).at(1)));
    }

    // return the evaluated contraction
    return std::stod(contractions.at(2)) * views.at(0).item<double>();
}

std::string MollerPlesset::get_name() const {
    return "MP" + std::to_string(input.order);
}

double MollerPlesset::run(const System& system, const torch::Tensor& F_MS, const torch::Tensor& J_MS_AP) const {
    // define the correlation energy container
    std::vector<double> energies(input.order - 1, 0);

    // loop over the orders of MP perturbation theory
    for (int i = 2; i <= input.order; i++) {

        // start the timer
        Timepoint order_timer = Timer::Now();

        // define contributions with the container for the evaluated contributions to the current order
        std::vector<std::string> contributions = Split(mbptf.at(i - 2), ' '); std::vector<double> Empi(contributions.size());

        // print the header
        std::printf("\nMP%02d ENERGY CALCULATION\n%13s %17s %12s\n", i, "CONTR", "VALUE", "TIME");

        // loop over the contractions
        for (size_t j = 0; j < contributions.size(); j++) {

            // start the timer
            Timepoint contribution_timer = Timer::Now();

            // add the contribution to the evaluated contribution container
            double Empii = evaluate_contraction(system, contributions.at(j), F_MS, J_MS_AP, i); Empi.at(j) = Empii;

            // print the contribution value with timing
            std::printf("%06d/%06d %17.14f %s\n", (int)j + 1, (int)contributions.size(), Empii, Timer::Format(Timer::Elapsed(contribution_timer)).c_str());
        }

        // print the MPN energy
        energies.at(i - 2) = std::accumulate(Empi.begin(), Empi.end(), 0.0); std::printf("MP%02d CORRELATION ENERGY: %.16f\n", i, energies.at(i - 2));

        // print the MPN timing
        std::printf("MP%02d CORRELATION ENERGY TIMING: %s\n", i, Timer::Format(Timer::Elapsed(order_timer)).c_str());
    }

    // return the correlation energy
    return std::accumulate(energies.begin(), energies.end(), 0.0);
}
