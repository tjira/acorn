#pragma once

#include     <fftw3.h>
#include <Eigen/Dense>

namespace FourierTransform {
    // general FFT functions
    void IFFT(std::complex<double>* in, std::complex<double>* out, const std::vector<int>& sizes);
    void FFT(std::complex<double>* in, std::complex<double>* out, const std::vector<int>& sizes);
    void C2RFFT(std::complex<double>* in, double* out, const std::vector<int>& sizes);

    // forward and inverse Fourier transforms
    Eigen::VectorXcd IFFT(const Eigen::VectorXcd& data, const std::vector<int>& sizes = {});
    Eigen::VectorXcd FFT(const Eigen::VectorXcd& data, const std::vector<int>& sizes = {});
    Eigen::VectorXd C2RFFT(const Eigen::VectorXcd& data, const std::vector<int>& sizes = {});
}
