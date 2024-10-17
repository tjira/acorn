#include "integral.h"

torch::Tensor Integral::double_electron(libint2::Engine& engine, const libint2::BasisSet& shells) {
    // define the number of basis functions, matrix of integrals and shell to basis function map
    int nbf = shells.nbf(); std::vector<double> ints(nbf * nbf * nbf * nbf, 0); std::vector<size_t> sh2bf = shells.shell2bf();

    // generate engine for every thread
    std::vector<libint2::Engine> engines(nthread, engine);

    // loop over all unique elements
    #pragma omp parallel for num_threads(nthread)
    for (size_t i = 0; i < shells.size(); i++) for (size_t j = i; j < shells.size(); j++) for (size_t k = i; k < shells.size(); k++) for (size_t l = (i == k ? j : k); l < shells.size(); l++) {

        // calculate the integral and skip if it is zero
        int integral_index = 0; engines.at(omp_get_thread_num()).compute(shells.at(i), shells.at(j), shells.at(k), shells.at(l)); if (engines.at(omp_get_thread_num()).results().at(0) == nullptr) continue;

        // assign the integrals
        for (size_t m = 0; m < shells.at(i).size(); m++) {
            for (size_t n = 0; n < shells.at(j).size(); n++) {
                for (size_t o = 0; o < shells.at(k).size(); o++) {
                    for (size_t p = 0; p < shells.at(l).size(); p++) {
                        ints.at((m + sh2bf.at(i)) * nbf * nbf * nbf + (n + sh2bf.at(j)) * nbf * nbf + (o + sh2bf.at(k)) * nbf + (p + sh2bf.at(l))) = engines.at(omp_get_thread_num()).results().at(0)[integral_index  ];
                        ints.at((m + sh2bf.at(i)) * nbf * nbf * nbf + (n + sh2bf.at(j)) * nbf * nbf + (p + sh2bf.at(l)) * nbf + (o + sh2bf.at(k))) = engines.at(omp_get_thread_num()).results().at(0)[integral_index  ];
                        ints.at((n + sh2bf.at(j)) * nbf * nbf * nbf + (m + sh2bf.at(i)) * nbf * nbf + (o + sh2bf.at(k)) * nbf + (p + sh2bf.at(l))) = engines.at(omp_get_thread_num()).results().at(0)[integral_index  ];
                        ints.at((n + sh2bf.at(j)) * nbf * nbf * nbf + (m + sh2bf.at(i)) * nbf * nbf + (p + sh2bf.at(l)) * nbf + (o + sh2bf.at(k))) = engines.at(omp_get_thread_num()).results().at(0)[integral_index  ];
                        ints.at((o + sh2bf.at(k)) * nbf * nbf * nbf + (p + sh2bf.at(l)) * nbf * nbf + (m + sh2bf.at(i)) * nbf + (n + sh2bf.at(j))) = engines.at(omp_get_thread_num()).results().at(0)[integral_index  ];
                        ints.at((o + sh2bf.at(k)) * nbf * nbf * nbf + (p + sh2bf.at(l)) * nbf * nbf + (n + sh2bf.at(j)) * nbf + (m + sh2bf.at(i))) = engines.at(omp_get_thread_num()).results().at(0)[integral_index  ];
                        ints.at((p + sh2bf.at(l)) * nbf * nbf * nbf + (o + sh2bf.at(k)) * nbf * nbf + (m + sh2bf.at(i)) * nbf + (n + sh2bf.at(j))) = engines.at(omp_get_thread_num()).results().at(0)[integral_index  ];
                        ints.at((p + sh2bf.at(l)) * nbf * nbf * nbf + (o + sh2bf.at(k)) * nbf * nbf + (n + sh2bf.at(j)) * nbf + (m + sh2bf.at(i))) = engines.at(omp_get_thread_num()).results().at(0)[integral_index++];
                    }
                }
            }
        }
    }

    // return the integrals
    return torch::from_blob(ints.data(), {nbf, nbf, nbf, nbf}, torch::kDouble).clone();
}

torch::Tensor Integral::single_electron(libint2::Engine& engine, const libint2::BasisSet& shells) {
    // define the number of basis functions, matrix of integrals and shell to basis function map
    int nbf = shells.nbf(); std::vector<double> ints(nbf * nbf, 0); std::vector<size_t> sh2bf = shells.shell2bf();

    // generate engine for every thread
    std::vector<libint2::Engine> engines(nthread, engine);

    // loop over all unique elements
    #pragma omp parallel for num_threads(nthread)
    for (size_t i = 0; i < shells.size(); i++) for (size_t j = i; j < shells.size(); j++) {

        // calculate the integral and skip if it is zero
        int integral_index = 0; engines.at(omp_get_thread_num()).compute(shells.at(i), shells.at(j)); if (engines.at(omp_get_thread_num()).results().at(0) == nullptr) continue;

        // assign the integrals
        for (size_t k = 0; k < shells.at(i).size(); k++) {
            for (size_t l = 0; l < shells.at(j).size(); l++) {
                ints.at((k + sh2bf.at(i)) * nbf + l + sh2bf.at(j)) = engines.at(omp_get_thread_num()).results().at(0)[integral_index  ];
                ints.at((l + sh2bf.at(j)) * nbf + k + sh2bf.at(i)) = engines.at(omp_get_thread_num()).results().at(0)[integral_index++];
            }
        }
    }

    // return the integrals
    return torch::from_blob(ints.data(), {nbf, nbf}, torch::kDouble).clone();
}

torch::Tensor Integral::double_electron_d1(libint2::Engine& engine, const libint2::BasisSet& shells, const std::vector<libint2::Atom>&) {
    // define the number of basis functions, matrix of integrals and shell to basis function map
    int nbf = shells.nbf(), dim = 6; std::vector<double> ints; std::vector<size_t> sh2bf = shells.shell2bf();

    // generate engine for every thread and initialize the integrals
    std::vector<libint2::Engine> engines(nthread, engine); ints.resize(nbf * nbf * nbf * nbf * dim, 0);

    // loop over all unique elements
    #pragma omp parallel for num_threads(nthread)
    for (size_t i = 0; i < shells.size(); i++) for (size_t j = 0; j < shells.size(); j++) for (size_t k = 0; k < shells.size(); k++) for (size_t l = 0; l < shells.size(); l++) {

        // calculate the integral and skip if it is zero
        int integral_index = 0; engines.at(omp_get_thread_num()).compute(shells.at(i), shells.at(j), shells.at(k), shells.at(l)); if (engines.at(omp_get_thread_num()).results().at(0) == nullptr) continue;

        // for every shell combination
        for (size_t m = 0; m < shells.at(i).size(); m++) {
            for (size_t n = 0; n < shells.at(j).size(); n++) {
                for (size_t o = 0; o < shells.at(k).size(); o++) {
                    for (size_t p = 0; p < shells.at(l).size(); p++) {

                        // for every dimension
                        for (int q = 0; q < dim; q++) {
                            ints.at((m + sh2bf.at(i)) * nbf * nbf * nbf * dim + (n + sh2bf.at(j)) * nbf * nbf * dim + (o + sh2bf.at(k)) * nbf * dim + (p + sh2bf.at(l)) * dim + q) = engines.at(omp_get_thread_num()).results().at(q)[integral_index];
                        }

                        // increment index
                        integral_index++;
                    }
                }
            }
        }
    }

    // return the integrals
    return torch::from_blob(ints.data(), {nbf, nbf, nbf, nbf, dim}, torch::kDouble).clone();
}

torch::Tensor Integral::single_electron_d1(libint2::Engine& engine, const libint2::BasisSet& shells, const std::vector<libint2::Atom>& atoms) {
    // define the number of basis functions and dimensions, matrix of integrals and shell to basis function map
    int nbf = shells.nbf(), dim = 6; std::vector<double> ints; std::vector<size_t> sh2bf = shells.shell2bf();

    // add dimensions if other coordinates are present
    if (engine.oper() == libint2::Operator::nuclear) dim += 3 * atoms.size();

    // generate engine for every thread and initialize the integrals
    std::vector<libint2::Engine> engines(nthread, engine); ints.resize(nbf * nbf * dim, 0);

    // loop over all unique elements
    #pragma omp parallel for num_threads(nthread)
    for (size_t i = 0; i < shells.size(); i++) for (size_t j = 0; j < shells.size(); j++) {

        // calculate the integral and skip if it is zero
        int integral_index = 0; engines.at(omp_get_thread_num()).compute(shells.at(i), shells.at(j)); if (engines.at(omp_get_thread_num()).results().at(0) == nullptr) continue;

        // for every shell pair
        for (size_t k = 0; k < shells.at(i).size(); k++) {
            for (size_t l = 0; l < shells.at(j).size(); l++) {

                // for every dimension
                for (int m = 0; m < dim; m++) {
                    ints.at((k + sh2bf.at(i)) * nbf * dim + (l + sh2bf.at(j)) * dim + m) = engines.at(omp_get_thread_num()).results().at(m)[integral_index];
                }

                // increment index
                integral_index++;
            }
        }
    }

    // return the integrals
    return torch::from_blob(ints.data(), {nbf, nbf, dim}, torch::kDouble).clone();
}

std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor> Integral::calculate(const std::vector<libint2::Atom>& atoms, const libint2::BasisSet& shells) const {
    // initialize libint
    libint2::initialize();

    // define the engines
    libint2::Engine kinetic_engine(libint2::Operator::kinetic, shells.max_nprim(), shells.max_l(), 0, precision);
    libint2::Engine nuclear_engine(libint2::Operator::nuclear, shells.max_nprim(), shells.max_l(), 0, precision);
    libint2::Engine overlap_engine(libint2::Operator::overlap, shells.max_nprim(), shells.max_l(), 0, precision);
    libint2::Engine coulomb_engine(libint2::Operator::coulomb, shells.max_nprim(), shells.max_l(), 0, precision);

    // set charges to the nuclear engine
    nuclear_engine.set_params(libint2::make_point_charges(atoms));
    
    // calculate the integrals
    torch::Tensor T = single_electron(kinetic_engine, shells);
    torch::Tensor V = single_electron(nuclear_engine, shells);
    torch::Tensor S = single_electron(overlap_engine, shells);
    torch::Tensor J = double_electron(coulomb_engine, shells);

    // finalize libint
    libint2::finalize();

    // return the integrals
    return std::make_tuple(T, V, S, J);
}

std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor> Integral::calculate_d1(const std::vector<libint2::Atom>& atoms, const libint2::BasisSet& shells) const {
    // initialize libint
    libint2::initialize();

    // define the engines
    libint2::Engine kinetic_engine_d1(libint2::Operator::kinetic, shells.max_nprim(), shells.max_l(), 1, precision);
    libint2::Engine nuclear_engine_d1(libint2::Operator::nuclear, shells.max_nprim(), shells.max_l(), 1, precision);
    libint2::Engine overlap_engine_d1(libint2::Operator::overlap, shells.max_nprim(), shells.max_l(), 1, precision);
    libint2::Engine coulomb_engine_d1(libint2::Operator::coulomb, shells.max_nprim(), shells.max_l(), 1, precision);

    // set charges to the nuclear engine
    nuclear_engine_d1.set_params(libint2::make_point_charges(atoms));
    
    // calculate the integrals
    torch::Tensor dT = single_electron_d1(kinetic_engine_d1, shells, atoms);
    torch::Tensor dV = single_electron_d1(nuclear_engine_d1, shells, atoms);
    torch::Tensor dS = single_electron_d1(overlap_engine_d1, shells, atoms);
    torch::Tensor dJ = double_electron_d1(coulomb_engine_d1, shells, atoms);

    // finalize libint
    libint2::finalize();

    // return the integrals
    return std::make_tuple(dT, dV, dS, dJ);
}
