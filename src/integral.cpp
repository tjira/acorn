#include "integral.h"

// include the libint
#include <libint2.hpp>

Matrix<> Integral::Single(libint2::Engine& engine, const System& system) {
    // extract shells, create a map between shell and basis function indices and initialize the engine vector
    auto shells = system.getShells(); std::vector<size_t> sh2bf = shells.shell2bf(); std::vector<libint2::Engine> engines(nthread, engine);

    // create the integral matrix
    Matrix<> matrix(shells.nbf(), shells.nbf());

    #pragma omp parallel for num_threads(nthread)
    for (size_t i = 0; i < shells.size(); i++) {
        for (size_t j = 0; j <= i; j++) {
            // calculate the integral
            engines.at(omp_get_thread_num()).compute(shells.at(j), shells.at(i));

            // obtain the results and skip if empty
            const auto& result = engines.at(omp_get_thread_num()).results(); if (result.at(0) == nullptr) continue;

            // extract the integral result to a buffer
            Eigen::Map<const Matrix<>> buffer(result.at(0), shells.at(i).size(), shells.at(j).size());

            // assign the buffer to the correct position in the result matrix
            matrix.block(sh2bf.at(i), sh2bf.at(j), shells.at(i).size(), shells.at(j).size()) = buffer;

            // if the element is not on the diagonal fill also the mirrored element
            if (i != j) {
                matrix.block(sh2bf.at(j), sh2bf.at(i), shells.at(j).size(), shells.at(i).size()) = buffer.transpose();
            }
        }
    }

    // return the integrals
    return matrix;
};

Tensor<> Integral::Double(libint2::Engine& engine, const System& system) {
    // extract system shells, create a map between shell and basis function indices and initialize the engine vector
    auto shells = system.getShells(); std::vector<size_t> sh2bf = shells.shell2bf(); std::vector<libint2::Engine> engines(nthread, engine);

    // create the result tensor and set the tensor elements to zero
    Tensor<> tensor(shells.nbf(), shells.nbf(), shells.nbf(), shells.nbf()); tensor.setZero();

    #pragma omp parallel for num_threads(nthread)
    for (size_t i = 0; i < shells.size(); i++) {
        for (size_t j = 0; j <= i; j++) {
            for (size_t k = 0; k <= i; k++) {
                for (size_t l = 0; l <= (i == k ? j : k); l++) {
                    // calculate the integral
                    engines.at(omp_get_thread_num()).compute(shells.at(i), shells.at(j), shells.at(k), shells.at(l));

                    // obtain the integral values and skip if empty
                    const auto& result = engines.at(omp_get_thread_num()).results(); if (result[0] == nullptr) continue;

                    // extract indices of the current basis functions
                    size_t bfi = sh2bf.at(i), bfj = sh2bf.at(j), bfk = sh2bf.at(k), bfl = sh2bf.at(l);

                    // assign the elements using the 8-fold symmetry
                    for (size_t m = 0, q = 0; m < shells.at(i).size(); m++) {
                        for (size_t n = 0; n < shells.at(j).size(); n++) {
                            for (size_t o = 0; o < shells.at(k).size(); o++) {
                                for (size_t p = 0; p < shells.at(l).size(); p++, q++) {
                                    tensor(m + bfi, n + bfj, o + bfk, p + bfl) = result.at(0)[q];
                                    tensor(m + bfi, n + bfj, p + bfl, o + bfk) = result.at(0)[q];
                                    tensor(n + bfj, m + bfi, o + bfk, p + bfl) = result.at(0)[q];
                                    tensor(n + bfj, m + bfi, p + bfl, o + bfk) = result.at(0)[q];
                                    tensor(o + bfk, p + bfl, m + bfi, n + bfj) = result.at(0)[q];
                                    tensor(o + bfk, p + bfl, n + bfj, m + bfi) = result.at(0)[q];
                                    tensor(p + bfl, o + bfk, m + bfi, n + bfj) = result.at(0)[q];
                                    tensor(p + bfl, o + bfk, n + bfj, m + bfi) = result.at(0)[q];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // return the integrals
    return tensor;
};

Tensor<> Integral::Coulomb(const System& system) {
    libint2::Engine engine(libint2::Operator::coulomb, system.getShells().max_nprim(), system.getShells().max_l(), 0, 1e-12);
    return Double(engine, system);
}

Matrix<> Integral::Kinetic(const System& system) {
    libint2::Engine engine(libint2::Operator::kinetic, system.getShells().max_nprim(), system.getShells().max_l(), 0, 1e-12);
    return Single(engine, system);
}

Matrix<> Integral::Nuclear(const System& system) {
    libint2::Engine engine(libint2::Operator::nuclear, system.getShells().max_nprim(), system.getShells().max_l(), 0, 1e-12);
    engine.set_params(libint2::make_point_charges(system.getAtoms<libint2::Atom>())); return Single(engine, system);
}

Matrix<> Integral::Overlap(const System& system) {
    libint2::Engine engine(libint2::Operator::overlap, system.getShells().max_nprim(), system.getShells().max_l(), 0, 1e-12);
    return Single(engine, system);
}

void Integral::Initialize() {libint2::initialize();}
void Integral::Finalize() {libint2::finalize();}
