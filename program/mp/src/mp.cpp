#include "mp.h"

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

#define SPLIT(S, D) [](const std::string& s, char d) {std::vector<std::string> r; std::stringstream ss(s); std::string l; while (getline(ss, l, d)) {r.push_back(l);} return r;}(S, D)

double Acorn::MBPT::evaluate(const std::string& contrstr, const torch::Tensor& Jmsa, const torch::Tensor& Ems, int nos, int o) {
    // define the contribution and split it
    std::vector<std::string> contraction = SPLIT(contrstr, ';');

    // replace the plus sign with a comma
    std::replace(contraction.at(0).begin(), contraction.at(0).end(), '+', ',');
    std::replace(contraction.at(1).begin(), contraction.at(1).end(), '+', ',');
    
    // define the tensor slices with axes, array of orbital energy slices and vector of ones for reshaping
    std::vector<torch::Tensor> views; std::vector<Slice> axes(4); std::vector<torch::Tensor> epsts; std::vector<int64_t> ones = {-1};

    // extract the contraction indices and split them
    std::vector<std::string> contrinds = SPLIT(contraction.at(0), ',');

    // add the coulomb contributions
    for (int k = 0; k < o; k++) {
        for (int l = 0; l < 4; l++) {
            if (std::string(OCC).find(contrinds.at(k).at(l)) != std::string::npos) axes.at(l) = Slice(None, nos);
            if (std::string(VRT).find(contrinds.at(k).at(l)) != std::string::npos) axes.at(l) = Slice(nos, None);
        } views.push_back(Jmsa.index({axes.at(0), axes.at(1), axes.at(2), axes.at(3)}));
    }

    // add the energy denominators
    for (int k = o; k < contrinds.size(); k++, epsts.clear(), ones = {-1}) {
        for (int l = 0; l < contrinds.at(k).size(); l++) {
            if (std::string(OCC).find(contrinds.at(k).at(l)) != std::string::npos) {
                epsts.push_back(+1*Ems.index({Slice(None, nos)}).reshape(at::IntArrayRef(ones))), ones.push_back(1);
            }
        }
        for (int l = 0; l < contrinds.at(k).size(); l++) {
            if (std::string(VRT).find(contrinds.at(k).at(l)) != std::string::npos) {
                epsts.push_back(-1*Ems.index({Slice(nos, None)}).reshape(at::IntArrayRef(ones))), ones.push_back(1);
            }
        }
        views.push_back(1 / (std::accumulate(epsts.begin() + 1, epsts.end(), epsts.at(0))));
    }

    // define the contraction path
    std::vector<std::tuple<std::string, std::array<int, 2>>> path;

    // split the nodes and add them to the path
    for (const std::string& node: SPLIT(contraction.at(1), ':')) {
        path.push_back({SPLIT(node, '/').at(1), {std::stoi(SPLIT(SPLIT(node, '/').at(0), '-').at(0)), std::stoi(SPLIT(SPLIT(node, '/').at(0), '-').at(1))}});
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
    return std::stod(contraction.at(2)) * views.at(0).item<double>();
}

void Acorn::MBPT::run(const Options& opt, std::vector<timepoint>& timers) {
    // start the timer for integral loading
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of integral loading
    std::cout << "SYSTEM AND INTEGRALS IN MS BASIS READING: " << std::flush;

    // load the integrals in MS basis
    torch::Tensor Jms = torch::ReadTensor("J_MS.mat").squeeze();
    torch::Tensor Vms = torch::ReadTensor("V_MS.mat").squeeze();
    torch::Tensor Tms = torch::ReadTensor("T_MS.mat").squeeze();
    torch::Tensor Ems = torch::ReadTensor("E_MS.mat").squeeze();
    torch::Tensor N   = torch::ReadTensor("N.mat"   ).squeeze();

    // print the time for integral loading
    std::cout << eltime(timers.at(1)) << std::endl;

    // extract the number of occupied and virtual spinorbitals
    int nos = 2 * N.index({0}).item<int>(); int nvs = Ems.sizes().at(0) - nos;

    // define the energy containers
    std::vector<double> Empn; double E = 0;

    // initialize the antisymmetrized Coulomb integrals in Physicists' notation and the Hamiltonian matrix in MS basis
    torch::Tensor Jmsa = (Jms - Jms.permute({0, 3, 2, 1})).permute({0, 2, 1, 3}); torch::Tensor Hms = Vms + Tms;

    // calculate the HF energy
    for (int i = 0; i < nos; i++) {
        E += Hms.index({i, i}).item<double>(); for (int j = 0; j < nos; j++) E += 0.5 * Jmsa.index({i, j, i, j}).item<double>();
    }

    // print new line
    std::cout << std::endl;

    // loop over the orders of MP perturbation theory
    for (int i = 2; i <= opt.order; i++) {

        // start the timer
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // define contributions with energy
        std::vector<std::string> contributions = SPLIT(mbptf.at(i - 2), ' '); std::vector<double> Empi(contributions.size());

        // print the required number of contributions
        std::printf("MP%02d ENERGY CALCULATION\n%13s %17s %12s\n", i, "CONTR", "VALUE", "TIME");

        // loop over the contractions
        #pragma omp parallel for num_threads(Acorn::MBPT::nthread)
        for (size_t j = 0; j < contributions.size(); j++) {

            // start the timer
            timers.at(2) = std::chrono::high_resolution_clock().now();

            // add the contribution to the energy
            double Empii = Acorn::MBPT::evaluate(contributions.at(j), Jmsa, Ems, nos, i); Empi.at(j) = Empii;

            // print the contribution
            std::printf("%06d/%06d %17.14f %s\n", (int)j + 1, (int)contributions.size(), Empii, eltime(timers.at(2)).c_str());
        }

        // print the MPN energy
        Empn.push_back(std::accumulate(Empi.begin(), Empi.end(), 0.0)); std::printf("MP%02d CORRELATION ENERGY: %.16f\n", i, Empn.back());

        // print the MPN timing
        std::printf("MP%02d CORRELATION ENERGY TIMING: %s\n", i, eltime(timers.at(1)).c_str());

        // print new line
        if (i < opt.order) std::cout << std::endl;
    }

    // sum the MPN energies
    E += std::accumulate(Empn.begin(), Empn.end(), 0.0);

    // print the final energy
    std::printf("\nFINAL SINGLE POINT ENERGY: %.13f\n\nTOTAL TIME: %s\n", E + N.index({1}).item<double>(), eltime(timers.at(0)).c_str());
}
