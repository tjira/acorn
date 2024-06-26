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

Vector FourierTransform::C2RFFT(const ComplexVector& data, const std::vector<int>& sizes) {
    // extract vector size
    int m = data.rows();

    // initialize the output vector and perform the C2RFFT
    std::vector<double> out((data.size() - 1) * 2); FourierTransform::C2RFFT(const_cast<ComplexVector&>(data).data(), out.data(), sizes.empty() ? std::vector<int>{m} : sizes);

    // transform the output to the correct format and return
    return Eigen::Map<Vector>(out.data(), m);
}

ComplexVector FourierTransform::IFFT(const ComplexVector& data, const std::vector<int>& sizes) {
    // extract vector size
    int m = data.rows();

    // initialize the output vector and perform the IFFT
    std::vector<std::complex<double>> out(data.size()); IFFT(const_cast<ComplexVector&>(data).data(), out.data(), sizes.empty() ? std::vector<int>{m} : sizes);

    // transform the output to the correct format and return
    return Eigen::Map<ComplexVector>(out.data(), m);
}

ComplexVector FourierTransform::FFT(const ComplexVector& data, const std::vector<int>& sizes) {
    // extract vector size
    int m = data.rows();

    // initialize the output vector and perform the FFT
    std::vector<std::complex<double>> out(data.size()); FFT(const_cast<ComplexVector&>(data).data(), out.data(), sizes.empty() ? std::vector<int>{m} : sizes);

    // transform the output to the correct format and return
    return Eigen::Map<ComplexVector>(out.data(), m) / out.size();
}
