#include "fourier.h"
#include <fftw3.h>

void FourierTransform::FFT(std::complex<double>* in, std::complex<double>* out, const std::vector<int>& sizes, int sign) {
    // create the plan for the FFT
    fftw_plan plan = fftw_plan_dft(sizes.size(), sizes.data(), reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), sign, FFTW_ESTIMATE);

    // execute the plan and destroy it
    fftw_execute(plan); fftw_destroy_plan(plan);
}

EigenMatrix<std::complex<double>> FourierTransform::Forward(EigenMatrix<std::complex<double>> data) {
    // extract dimension
    int m = data.rows(), n = data.cols();

    // initialize the output vector and perform the FFT
    std::vector<std::complex<double>> out(data.size()); FFT(data.data(), out.data(), std::vector<int>{m, n}, 1);

    // transform the output to the correct format and return
    return Eigen::Map<EigenMatrix<std::complex<double>>>(out.data(), m, n) / out.size();
}

EigenMatrix<std::complex<double>> FourierTransform::Inverse(EigenMatrix<std::complex<double>> data) {
    // extract the dimensions
    int m = data.rows(), n = data.cols();

    // initialize the output vector and perform the FFT
    std::vector<std::complex<double>> out(data.size()); FFT(data.data(), out.data(), std::vector<int>{m, n}, -1);

    // transform the output to the correct format and return
    return Eigen::Map<EigenMatrix<std::complex<double>>>(out.data(), m, n);
}
