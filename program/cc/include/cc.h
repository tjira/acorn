#pragma once

#include "tensor.h"

using namespace torch::indexing;

namespace Acorn {
    namespace CC {
        namespace LCCD {
            torch::Tensor amplitude(const torch::Tensor& Jmsa, const torch::Tensor& Emsd, const torch::Tensor& T2, int nos);
            double energy(const torch::Tensor& Jmsa, const torch::Tensor& T2, int nos);
        }
        namespace CCD {
            torch::Tensor amplitude(const torch::Tensor& Jmsa, const torch::Tensor& Emsd, const torch::Tensor& T2, int nos);
            double energy(const torch::Tensor& Jmsa, const torch::Tensor& T2, int nos);
        }
        namespace CCSD {
            std::tuple<torch::Tensor, torch::Tensor> amplitude(const torch::Tensor& Jmsa, const torch::Tensor& Fms, const torch::Tensor& Emss, const torch::Tensor& Emsd, const torch::Tensor& T1, const torch::Tensor& T2, int nos);
            double pertrubationTriple(const torch::Tensor& Jmsa, const torch::Tensor& Emst, const torch::Tensor& T1, const torch::Tensor& T2, int nos);
            double energy(const torch::Tensor& Jmsa, const torch::Tensor& Fms, const torch::Tensor& T1, const torch::Tensor& T2, int nos);
        }
    }
}
