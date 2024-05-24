#pragma once

#include "linalg.h"

namespace FourierTransform {
    // general FFT function
    void FFT(std::complex<double>* in, std::complex<double>* out, const std::vector<int>& sizes, int sign);

    // forward and inverse Fourier transforms
    EigenMatrix<std::complex<double>> Forward(EigenMatrix<std::complex<double>> data);
    EigenMatrix<std::complex<double>> Inverse(EigenMatrix<std::complex<double>> data);
}
