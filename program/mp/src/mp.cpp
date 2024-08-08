#include "mp.h"

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
