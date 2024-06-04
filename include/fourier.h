#pragma once

#include "linalg.h"

namespace FourierTransform {
    // general FFT functions
    void IFFT(std::complex<double>* in, std::complex<double>* out, const std::vector<int>& sizes);
    void FFT(std::complex<double>* in, std::complex<double>* out, const std::vector<int>& sizes);
    void C2RFFT(std::complex<double>* in, double* out, const std::vector<int>& sizes);

    // forward and inverse Fourier transforms
    ComplexMatrix IFFT(const ComplexMatrix& data);
    ComplexMatrix FFT(const ComplexMatrix& data);
    Matrix C2RFFT(const ComplexMatrix& data);
}
