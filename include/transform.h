#pragma once

#include <torch/torch.h>

namespace Transform {
    // two electron integrals
    torch::Tensor DoubleSpin(const torch::Tensor& J_AO, const torch::Tensor& C_MO);
    torch::Tensor DoubleSpatial(const torch::Tensor& J_AO, torch::Tensor& C_MO);

    // antisymmetric two electron integrals in physicist notation
    torch::Tensor DoubleSpinAntsymPhys(const torch::Tensor& J_AO, const torch::Tensor& C_MO);
    torch::Tensor DoubleSpatialAntsymPhys(const torch::Tensor& J_AO, torch::Tensor& C_MO);

    // excitation energy fraction used in post-HF methods
    torch::Tensor ExcitationEnergyFraction(const torch::Tensor& F_MS, at::indexing::Slice o, at::indexing::Slice v, int order);

    // one electron integrals
    torch::Tensor SingleSpin(const torch::Tensor& A_AO, const torch::Tensor& C_MO);
    torch::Tensor SingleSpatial(const torch::Tensor& A_AO, torch::Tensor& C_MO);
}
