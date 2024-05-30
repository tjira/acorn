#include "fourier.h"
#include <fftw3.h>

void FourierTransform::C2RFFT(std::complex<double>* in, double* out, const std::vector<int>& sizes) {
    // create the plan for the C2RFFT
    fftw_plan plan = fftw_plan_dft_c2r(sizes.size(), sizes.data(), reinterpret_cast<fftw_complex*>(in), out, FFTW_ESTIMATE);

    // execute the plan and destroy it
    fftw_execute(plan); fftw_destroy_plan(plan);
}

void FourierTransform::IFFT(std::complex<double>* in, std::complex<double>* out, const std::vector<int>& sizes) {
    // create the plan for the IFFT
    fftw_plan plan = fftw_plan_dft(sizes.size(), sizes.data(), reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), -1, FFTW_ESTIMATE);

    // execute the plan and destroy it
    fftw_execute(plan); fftw_destroy_plan(plan);
}

void FourierTransform::FFT(std::complex<double>* in, std::complex<double>* out, const std::vector<int>& sizes) {
    // create the plan for the FFT
    fftw_plan plan = fftw_plan_dft(sizes.size(), sizes.data(), reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), 1, FFTW_ESTIMATE);

    // execute the plan and destroy it
    fftw_execute(plan); fftw_destroy_plan(plan);
}

EigenMatrix<double> FourierTransform::C2RFFT(const EigenMatrix<std::complex<double>>& data) {
    // extract dimension
    int m = data.rows(), n = data.cols();

    // initialize the output vector and perform the C2RFFT
    std::vector<double> out((data.size() - 1) * 2); FourierTransform::C2RFFT(const_cast<EigenMatrix<std::complex<double>>&>(data).data(), out.data(), std::vector<int>{m, n});

    // transform the output to the correct format and return
    return Eigen::Map<EigenVector<double>>(out.data(), data.size());
}

EigenMatrix<std::complex<double>> FourierTransform::IFFT(const EigenMatrix<std::complex<double>>& data) {
    // extract the dimensions
    int m = data.rows(), n = data.cols();

    // initialize the output vector and perform the IFFT
    std::vector<std::complex<double>> out(data.size()); IFFT(const_cast<EigenMatrix<std::complex<double>>&>(data).data(), out.data(), std::vector<int>{m, n});

    // transform the output to the correct format and return
    return Eigen::Map<EigenMatrix<std::complex<double>>>(out.data(), m, n);
}

EigenMatrix<std::complex<double>> FourierTransform::FFT(const EigenMatrix<std::complex<double>>& data) {
    // extract dimension
    int m = data.rows(), n = data.cols();

    // initialize the output vector and perform the FFT
    std::vector<std::complex<double>> out(data.size()); FFT(const_cast<EigenMatrix<std::complex<double>>&>(data).data(), out.data(), std::vector<int>{m, n});

    // transform the output to the correct format and return
    return Eigen::Map<EigenMatrix<std::complex<double>>>(out.data(), m, n) / out.size();
}
