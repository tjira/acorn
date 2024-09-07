#pragma once

#include <torch/torch.h>

namespace Transform {
    // MO to MS transformation for coefficient matrix
    torch::Tensor CoefficientSpin(const torch::Tensor& C_MO);

    // two electron integrals
    torch::Tensor DoubleSpatial(const torch::Tensor& J_AO, const torch::Tensor& C_MO);
    torch::Tensor DoubleSpin(const torch::Tensor& J_AO, const torch::Tensor& C_MS);

    // antisymmetric two electron integrals in physicist notation
    torch::Tensor DoubleSpatialAntsymPhys(const torch::Tensor& J_AO, const torch::Tensor& C_MO);
    torch::Tensor DoubleSpinAntsymPhys(const torch::Tensor& J_AO, const torch::Tensor& C_MS);

    // excitation energy fraction used in post-HF methods
    torch::Tensor ExcitationEnergyFraction(const torch::Tensor& F_MS, at::indexing::Slice o, at::indexing::Slice v, int order);

    // one electron integrals
    torch::Tensor SingleSpatial(const torch::Tensor& A_AO, const torch::Tensor& C_MO);
    torch::Tensor SingleSpin(const torch::Tensor& A_AO, const torch::Tensor& C_MS);
}
