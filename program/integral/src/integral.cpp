#include "integral.h"

torch::Tensor Integral::Single(libint2::Engine& engine, const libint2::BasisSet& shells) {
    // define the number of basis functions, matrix of integrals and shell to basis function map
    int nbf = shells.nbf(); std::vector<double> ints(nbf * nbf, 0); std::vector<size_t> sh2bf = shells.shell2bf();

    // loop over all unique elements
    for (size_t i = 0; i < shells.size(); i++) {
        for (size_t j = i; j < shells.size(); j++) {

            // calculate the integral and skip if it is zero
            engine.compute(shells.at(i), shells.at(j)); if (engine.results().at(0) == nullptr) continue;

            // assign the integrals
            for (size_t k = 0, m = 0; k < shells.at(i).size(); k++) {
                for (size_t l = 0; l < shells.at(j).size(); l++, m++) {
                    ints.at((k + sh2bf.at(i)) * nbf + l + sh2bf.at(j)) = engine.results().at(0)[m];
                    ints.at((l + sh2bf.at(j)) * nbf + k + sh2bf.at(i)) = engine.results().at(0)[m];
                }
            }
        }
    }

    // return the integrals
    return torch::from_blob(ints.data(), {nbf, nbf}, torch::kDouble).clone();
}

torch::Tensor Integral::Double(libint2::Engine& engine, const libint2::BasisSet& shells) {
    // define the number of basis functions, matrix of integrals and shell to basis function map
    int nbf = shells.nbf(); std::vector<double> ints(nbf * nbf * nbf * nbf, 0); std::vector<size_t> sh2bf = shells.shell2bf();

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
                                    ints.at((m + sh2bf.at(i)) * nbf * nbf * nbf + (n + sh2bf.at(j)) * nbf * nbf + (o + sh2bf.at(k)) * nbf + (p + sh2bf.at(l))) = engine.results().at(0)[q];
                                    ints.at((m + sh2bf.at(i)) * nbf * nbf * nbf + (n + sh2bf.at(j)) * nbf * nbf + (p + sh2bf.at(l)) * nbf + (o + sh2bf.at(k))) = engine.results().at(0)[q];
                                    ints.at((n + sh2bf.at(j)) * nbf * nbf * nbf + (m + sh2bf.at(i)) * nbf * nbf + (o + sh2bf.at(k)) * nbf + (p + sh2bf.at(l))) = engine.results().at(0)[q];
                                    ints.at((n + sh2bf.at(j)) * nbf * nbf * nbf + (m + sh2bf.at(i)) * nbf * nbf + (p + sh2bf.at(l)) * nbf + (o + sh2bf.at(k))) = engine.results().at(0)[q];
                                    ints.at((o + sh2bf.at(k)) * nbf * nbf * nbf + (p + sh2bf.at(l)) * nbf * nbf + (m + sh2bf.at(i)) * nbf + (n + sh2bf.at(j))) = engine.results().at(0)[q];
                                    ints.at((o + sh2bf.at(k)) * nbf * nbf * nbf + (p + sh2bf.at(l)) * nbf * nbf + (n + sh2bf.at(j)) * nbf + (m + sh2bf.at(i))) = engine.results().at(0)[q];
                                    ints.at((p + sh2bf.at(l)) * nbf * nbf * nbf + (o + sh2bf.at(k)) * nbf * nbf + (m + sh2bf.at(i)) * nbf + (n + sh2bf.at(j))) = engine.results().at(0)[q];
                                    ints.at((p + sh2bf.at(l)) * nbf * nbf * nbf + (o + sh2bf.at(k)) * nbf * nbf + (n + sh2bf.at(j)) * nbf + (m + sh2bf.at(i))) = engine.results().at(0)[q];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // return the integrals
    return torch::from_blob(ints.data(), {nbf, nbf, nbf, nbf}, torch::kDouble).clone();
}

torch::Tensor Integral::Nuclear(const std::vector<libint2::Atom>& atoms, const libint2::BasisSet& shells) {
    libint2::Engine engine(libint2::Operator::nuclear, shells.max_nprim(), shells.max_l(), 0, 1e-12);
    engine.set_params(libint2::make_point_charges(atoms)); return Single(engine, shells);
}

torch::Tensor Integral::Kinetic(const libint2::BasisSet& shells) {
    libint2::Engine engine(libint2::Operator::kinetic, shells.max_nprim(), shells.max_l(), 0, 1e-12);
    return Single(engine, shells);
}

torch::Tensor Integral::Overlap(const libint2::BasisSet& shells) {
    libint2::Engine engine(libint2::Operator::overlap, shells.max_nprim(), shells.max_l(), 0, 1e-12);
    return Single(engine, shells);
}

torch::Tensor Integral::Coulomb(const libint2::BasisSet& shells) {
    libint2::Engine engine(libint2::Operator::coulomb, shells.max_nprim(), shells.max_l(), 0, 1e-12);
    return Double(engine, shells);
}
