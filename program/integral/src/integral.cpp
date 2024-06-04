#include "integral.h"

Matrix Integral::Single(libint2::Engine& engine, const libint2::BasisSet& shells) {
    // define the number of basis functions, matrix of integrals and shell to basis function map
    int nbf = shells.nbf(); Matrix ints(nbf, nbf); ints.setZero(); std::vector<size_t> sh2bf = shells.shell2bf();

    // loop over all unique elements
    for (size_t i = 0; i < shells.size(); i++) {
        for (size_t j = i; j < shells.size(); j++) {

            // calculate the integral and skip if it is zero
            engine.compute(shells.at(i), shells.at(j)); if (engine.results().at(0) == nullptr) continue;

            // assign the integrals
            for (size_t k = 0, m = 0; k < shells.at(i).size(); k++) {
                for (size_t l = 0; l < shells.at(j).size(); l++, m++) {
                    ints(k + sh2bf.at(i), l + sh2bf.at(j)) = engine.results().at(0)[m];
                    ints(l + sh2bf.at(j), k + sh2bf.at(i)) = engine.results().at(0)[m];
                }
            }
        }
    }

    // return the integrals
    return ints;
}

Tensor<4> Integral::Double(libint2::Engine& engine, const libint2::BasisSet& shells) {
    // define the number of basis functions, matrix of integrals and shell to basis function map
    int nbf = shells.nbf(); Tensor<4> ints(nbf, nbf, nbf, nbf); ints.setZero(); std::vector<size_t> sh2bf = shells.shell2bf();

    // loop over all unique elements
    for (size_t i = 0; i < shells.size(); i++) {
        for (size_t j = i; j < shells.size(); j++) {
            for (size_t k = i; k < shells.size(); k++) {
                for (size_t l = (i == k ? j : k); l < shells.size(); l++) {

                    // calculate the integral, skip if it is zero
                    engine.compute(shells.at(i), shells.at(j), shells.at(k), shells.at(l)); if (engine.results().at(0) == nullptr) continue;

                    // assign the integrals
                    for (size_t m = 0, q = 0; m < shells.at(i).size(); m++) {
                        for (size_t n = 0; n < shells.at(j).size(); n++) {
                            for (size_t o = 0; o < shells.at(k).size(); o++) {
                                for (size_t p = 0; p < shells.at(l).size(); p++, q++) {
                                    ints(m + sh2bf.at(i), n + sh2bf.at(j), o + sh2bf.at(k), p + sh2bf.at(l)) = engine.results().at(0)[q];
                                    ints(m + sh2bf.at(i), n + sh2bf.at(j), p + sh2bf.at(l), o + sh2bf.at(k)) = engine.results().at(0)[q];
                                    ints(n + sh2bf.at(j), m + sh2bf.at(i), o + sh2bf.at(k), p + sh2bf.at(l)) = engine.results().at(0)[q];
                                    ints(n + sh2bf.at(j), m + sh2bf.at(i), p + sh2bf.at(l), o + sh2bf.at(k)) = engine.results().at(0)[q];
                                    ints(o + sh2bf.at(k), p + sh2bf.at(l), m + sh2bf.at(i), n + sh2bf.at(j)) = engine.results().at(0)[q];
                                    ints(o + sh2bf.at(k), p + sh2bf.at(l), n + sh2bf.at(j), m + sh2bf.at(i)) = engine.results().at(0)[q];
                                    ints(p + sh2bf.at(l), o + sh2bf.at(k), m + sh2bf.at(i), n + sh2bf.at(j)) = engine.results().at(0)[q];
                                    ints(p + sh2bf.at(l), o + sh2bf.at(k), n + sh2bf.at(j), m + sh2bf.at(i)) = engine.results().at(0)[q];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // return the integrals
    return ints;
}

Matrix Integral::Nuclear(const std::vector<libint2::Atom>& atoms, const libint2::BasisSet& shells) {
    libint2::Engine engine(libint2::Operator::nuclear, shells.max_nprim(), shells.max_l(), 0, 1e-12);
    engine.set_params(libint2::make_point_charges(atoms)); return Single(engine, shells);
}

Matrix Integral::Kinetic(const libint2::BasisSet& shells) {
    libint2::Engine engine(libint2::Operator::kinetic, shells.max_nprim(), shells.max_l(), 0, 1e-12);
    return Single(engine, shells);
}

Matrix Integral::Overlap(const libint2::BasisSet& shells) {
    libint2::Engine engine(libint2::Operator::overlap, shells.max_nprim(), shells.max_l(), 0, 1e-12);
    return Single(engine, shells);
}

Tensor<4> Integral::Coulomb(const libint2::BasisSet& shells) {
    libint2::Engine engine(libint2::Operator::coulomb, shells.max_nprim(), shells.max_l(), 0, 1e-12);
    return Double(engine, shells);
}
