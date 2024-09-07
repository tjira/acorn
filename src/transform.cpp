#include "transform.h"

torch::Tensor Transform::CoefficientSpin(const torch::Tensor& C_MO) {
    // create the tile matrix P
    torch::Tensor P = torch::eye(C_MO.sizes().at(0), C_MO.sizes().at(1), torch::kDouble).repeat_interleave(2, 1);

    // create the spin mask matrices
    torch::Tensor N = torch::arange(2, torch::kDouble).repeat({C_MO.sizes().at(0), C_MO.sizes().at(1)}); torch::Tensor M = (N + 1) % 2;

    // transform the wfn coefficients to the spin basis
    return torch::cat({C_MO.mm(P), C_MO.mm(P)}) * torch::cat({N, M});
}

torch::Tensor Transform::DoubleSpatial(const torch::Tensor& J, const torch::Tensor& C) {
    return torch::einsum("ip,jq,ijkl,kr,ls->pqrs", {C, C, J, C, C});
}

torch::Tensor Transform::DoubleSpin(const torch::Tensor& J_AO, const torch::Tensor& C_MS) {
    return DoubleSpatial(torch::kron(torch::eye(2, 2, torch::kDouble), torch::kron(torch::eye(2, 2, torch::kDouble), J_AO).swapaxes(0, 3).swapaxes(1, 2).contiguous()), C_MS);
}

torch::Tensor Transform::DoubleSpatialAntsymPhys(const torch::Tensor& J_AO, const torch::Tensor& C_MO) {
    torch::Tensor J_MO = DoubleSpatial(J_AO, C_MO); return (J_MO - J_MO.swapaxes(1, 3)).permute({0, 2, 1, 3});
}

torch::Tensor Transform::DoubleSpinAntsymPhys(const torch::Tensor& J_AO, const torch::Tensor& C_MS) {
    torch::Tensor J_MO = DoubleSpin(J_AO, C_MS); return (J_MO - J_MO.swapaxes(1, 3)).permute({0, 2, 1, 3});
}

torch::Tensor Transform::ExcitationEnergyFraction(const torch::Tensor& F_MS, at::indexing::Slice o, at::indexing::Slice v, int order) {
    // create the shaping vector and initialize the terms vector
    std::vector<int64_t> shape = {-1}; std::vector<torch::Tensor> terms;

    // add the terms to the tensor
    for (int i = 0; i < order; i++) terms.push_back(+1*F_MS.diagonal().index({o}).reshape(shape)), shape.push_back(1);
    for (int i = 0; i < order; i++) terms.push_back(-1*F_MS.diagonal().index({v}).reshape(shape)), shape.push_back(1);

    // return the tensor
    return 1 / std::accumulate(terms.begin() + 1, terms.end(), terms.at(0));
}

torch::Tensor Transform::SingleSpatial(const torch::Tensor& A_AO, const torch::Tensor& C_MO) {
    return C_MO.t().mm(A_AO).mm(C_MO);
}

torch::Tensor Transform::SingleSpin(const torch::Tensor& A_AO, const torch::Tensor& C_MS) {
    return SingleSpatial(torch::kron(torch::eye(2, 2, torch::kDouble), A_AO), C_MS);
}
