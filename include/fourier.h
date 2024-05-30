#pragma once

#include "linalg.h"

namespace FourierTransform {
    // general FFT functions
    void IFFT(std::complex<double>* in, std::complex<double>* out, const std::vector<int>& sizes);
    void FFT(std::complex<double>* in, std::complex<double>* out, const std::vector<int>& sizes);
    void C2RFFT(std::complex<double>* in, double* out, const std::vector<int>& sizes);

    // forward and inverse Fourier transforms
    EigenMatrix<std::complex<double>> IFFT(const EigenMatrix<std::complex<double>>& data);
    EigenMatrix<std::complex<double>> FFT(const EigenMatrix<std::complex<double>>& data);
    EigenMatrix<double> C2RFFT(const EigenMatrix<std::complex<double>>& data);
}
