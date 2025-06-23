#include <libint2.hpp>

std::pair<std::vector<libint2::Atom>, std::vector<libint2::Shell>> load(int natoms, const double *anums, const double *coords, int nbasis, const double *basis) {
    std::vector<libint2::Atom> atoms; std::vector<libint2::Shell> shells;

    for (int i = 0; i < natoms; i++) {
        atoms.push_back({(int)anums[i], coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]});
    }

    for (int i = 0; i < nbasis; i += 2 * basis[i] + 5) {

        int n = basis[i]; libint2::svector<double> alpha(n), c(n);

        for (int j = 0; j < basis[i]; j++) {
            alpha[j] = basis[i + 5 + j]; c[j] = basis[i + 5 + n + j];
        }

        shells.push_back(libint2::Shell{alpha, {{(int)basis[i + 1], false, c}}, {{basis[i + 2], basis[i + 3], basis[i + 4]}}});
    }

    return {atoms, shells};
}

#include <omp.h>
extern "C" {
    using namespace libint2; typedef unsigned long ulong;

    void oneelec(double *I, libint2::Engine &engine, const std::vector<Shell> &shells) {
        int nbf = 0; std::vector<size_t> sh2bf; std::vector<Engine> engines(std::getenv("OMP_NUM_THREADS") ? std::max(1, std::atoi(std::getenv("OMP_NUM_THREADS"))) : 1, engine);

        for (int i = 0, j = 0; j < shells.size(); i += shells[j].size(), j++) sh2bf.push_back(i);

        for (const auto& shell : shells) nbf += shell.size();

        #pragma omp parallel for num_threads(engines.size()) schedule(dynamic)
        for (size_t i = 0; i < shells.size(); i++) for (size_t j = i; j < shells.size(); j++) {

            int id = omp_get_thread_num(); int integral_index = 0; engines.at(id).compute(shells.at(i), shells.at(j)); if (engines.at(id).results().at(0) == nullptr) continue;

            for (size_t k = 0; k < shells.at(i).size(); k++) {
                for (size_t l = 0; l < shells.at(j).size(); l++) {
                    I[(k + sh2bf.at(i)) * nbf + (l + sh2bf.at(j))] = engines.at(id).results().at(0)[integral_index  ];
                    I[(l + sh2bf.at(j)) * nbf + (k + sh2bf.at(i))] = engines.at(id).results().at(0)[integral_index++];
                }
            }
        }
    }

    void twoelec(double *I, libint2::Engine &engine, const std::vector<Shell> &shells) {
        int nbf = 0; std::vector<size_t> sh2bf; std::vector<Engine> engines(std::getenv("OMP_NUM_THREADS") ? std::max(1, std::atoi(std::getenv("OMP_NUM_THREADS"))) : 1, engine);

        for (int i = 0, j = 0; j < shells.size(); i += shells[j].size(), j++) sh2bf.push_back(i);

        for (const auto& shell : shells) nbf += shell.size();

        #pragma omp parallel for num_threads(engines.size()) schedule(dynamic)
        for (size_t i = 0; i < shells.size(); i++) for (size_t j = i; j < shells.size(); j++) for (size_t k = i; k < shells.size(); k++) for (size_t l = (i == k ? j : k); l < shells.size(); l++) {

            int id = omp_get_thread_num(); int integral_index = 0; engines.at(id).compute(shells.at(i), shells.at(j), shells.at(k), shells.at(l)); if (engines.at(id).results().at(0) == nullptr) continue;

            for (size_t m = 0; m < shells.at(i).size(); m++) {
                for (size_t n = 0; n < shells.at(j).size(); n++) {
                    for (size_t o = 0; o < shells.at(k).size(); o++) {
                        for (size_t p = 0; p < shells.at(l).size(); p++) {

                            int bf1 = m + sh2bf.at(i); int bf2 = n + sh2bf.at(j); int bf3 = o + sh2bf.at(k); int bf4 = p + sh2bf.at(l);

                            I[bf1 * nbf * nbf * nbf + bf2 * nbf * nbf + bf3 * nbf + bf4] = engines.at(id).results().at(0)[integral_index  ];
                            I[bf1 * nbf * nbf * nbf + bf2 * nbf * nbf + bf4 * nbf + bf3] = engines.at(id).results().at(0)[integral_index  ];
                            I[bf2 * nbf * nbf * nbf + bf1 * nbf * nbf + bf3 * nbf + bf4] = engines.at(id).results().at(0)[integral_index  ];
                            I[bf2 * nbf * nbf * nbf + bf1 * nbf * nbf + bf4 * nbf + bf3] = engines.at(id).results().at(0)[integral_index  ];
                            I[bf3 * nbf * nbf * nbf + bf4 * nbf * nbf + bf1 * nbf + bf2] = engines.at(id).results().at(0)[integral_index  ];
                            I[bf3 * nbf * nbf * nbf + bf4 * nbf * nbf + bf2 * nbf + bf1] = engines.at(id).results().at(0)[integral_index  ];
                            I[bf4 * nbf * nbf * nbf + bf3 * nbf * nbf + bf1 * nbf + bf2] = engines.at(id).results().at(0)[integral_index  ];
                            I[bf4 * nbf * nbf * nbf + bf3 * nbf * nbf + bf2 * nbf + bf1] = engines.at(id).results().at(0)[integral_index++];
                        }
                    }
                }
            }
        }
    }

    void twoelec2fock(double *F, libint2::Engine &engine, const std::vector<Shell> &shells, const double *D) {
        int nbf = 0; std::vector<size_t> sh2bf; std::vector<Engine> engines(std::getenv("OMP_NUM_THREADS") ? std::max(1, std::atoi(std::getenv("OMP_NUM_THREADS"))) : 1, engine);

        for (int i = 0, j = 0; j < shells.size(); i += shells[j].size(), j++) sh2bf.push_back(i);

        for (const auto& shell : shells) nbf += shell.size();

        for (size_t i = 0; i < shells.size(); i++) for (size_t j = i; j < shells.size(); j++) for (size_t k = i; k < shells.size(); k++) for (size_t l = (i == k ? j : k); l < shells.size(); l++) {

            int id = omp_get_thread_num(); int integral_index = 0; engines.at(id).compute(shells.at(i), shells.at(j), shells.at(k), shells.at(l)); if (engines.at(id).results().at(0) == nullptr) continue;

            int degeneracy = (i==j ? 1 : 2) * (k==l ? 1 : 2) * (i==k ? (j==l ? 1 : 2) : 2);

            for (size_t m = 0; m < shells.at(i).size(); m++) {
                for (size_t n = 0; n < shells.at(j).size(); n++) {
                    for (size_t o = 0; o < shells.at(k).size(); o++) {
                        for (size_t p = 0; p < shells.at(l).size(); p++) {

                            int bf1 = m + sh2bf.at(i); int bf2 = n + sh2bf.at(j); int bf3 = o + sh2bf.at(k); int bf4 = p + sh2bf.at(l);

                            double value = engines.at(id).results().at(0)[integral_index++];

                            F[bf1 * nbf + bf2] += 0.500 * D[bf3 * nbf + bf4] * value * degeneracy;
                            F[bf3 * nbf + bf4] += 0.500 * D[bf1 * nbf + bf2] * value * degeneracy;
                            F[bf1 * nbf + bf3] -= 0.125 * D[bf2 * nbf + bf4] * value * degeneracy;
                            F[bf2 * nbf + bf4] -= 0.125 * D[bf1 * nbf + bf3] * value * degeneracy;
                            F[bf1 * nbf + bf4] -= 0.125 * D[bf2 * nbf + bf3] * value * degeneracy;
                            F[bf2 * nbf + bf3] -= 0.125 * D[bf1 * nbf + bf4] * value * degeneracy;
                        }
                    }
                }
            }
        }

        for (size_t i = 0; i < nbf; i++) for (size_t j = i + 1; j < nbf; j++) {
            F[i * nbf + j] = F[j * nbf + i] = (F[i * nbf + j] + F[j * nbf + i]) / 2;
        }
    }

    void twoelec2fockg(double *F, libint2::Engine &engine, const std::vector<Shell> &shells, const double *D) {
        int nbf = 0; std::vector<size_t> sh2bf; std::vector<Engine> engines(std::getenv("OMP_NUM_THREADS") ? std::max(1, std::atoi(std::getenv("OMP_NUM_THREADS"))) : 1, engine);

        for (int i = 0, j = 0; j < shells.size(); i += shells[j].size(), j++) sh2bf.push_back(i);

        for (const auto& shell : shells) nbf += 2 * shell.size();

        for (size_t i = 0; i < shells.size(); i++) for (size_t j = i; j < shells.size(); j++) for (size_t k = i; k < shells.size(); k++) for (size_t l = (i == k ? j : k); l < shells.size(); l++) {

            int id = omp_get_thread_num(); int integral_index = 0; engines.at(id).compute(shells.at(i), shells.at(j), shells.at(k), shells.at(l)); if (engines.at(id).results().at(0) == nullptr) continue;

            int degeneracy = (i == j ? 1 : 2) * (k == l ? 1 : 2) * (i == k ? (j == l ? 1 : 2) : 2);

            for (size_t m = 0; m < shells.at(i).size(); m++) {
                for (size_t n = 0; n < shells.at(j).size(); n++) {
                    for (size_t o = 0; o < shells.at(k).size(); o++) {
                        for (size_t p = 0; p < shells.at(l).size(); p++) {

                            int bf1 = m + sh2bf.at(i); int bf2 = n + sh2bf.at(j); int bf3 = o + sh2bf.at(k); int bf4 = p + sh2bf.at(l);

                            double value = engines.at(id).results().at(0)[integral_index++];

                            for (int s = 0; s < 2; s++) {
                                for (int t = 0; t < 2; t++) {

                                    int sbf1 = s * nbf / 2 + bf1; int sbf2 = t * nbf / 2 + bf2; int sbf3 = s * nbf / 2 + bf3; int sbf4 = t * nbf / 2 + bf4;

                                    F[sbf1 * nbf + sbf2] += D[sbf3 * nbf + sbf4] * value * degeneracy;
                                    F[sbf3 * nbf + sbf4] += D[sbf1 * nbf + sbf2] * value * degeneracy;

                                    if (s == t) {
                                        F[sbf1 * nbf + sbf3] -= 0.25 * D[sbf2 * nbf + sbf4] * value * degeneracy;
                                        F[sbf2 * nbf + sbf4] -= 0.25 * D[sbf1 * nbf + sbf3] * value * degeneracy;
                                        F[sbf1 * nbf + sbf4] -= 0.25 * D[sbf2 * nbf + sbf3] * value * degeneracy;
                                        F[sbf2 * nbf + sbf3] -= 0.25 * D[sbf1 * nbf + sbf4] * value * degeneracy;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (size_t i = 0; i < nbf; i++) for (size_t j = i + 1; j < nbf; j++) {
            F[i * nbf + j] = F[j * nbf + i] = (F[i * nbf + j] + F[j * nbf + i]) / 2;
        }
    }

    void coulomb(double *I, ulong natoms, const double *anums, const double *coords, int nbasis, const double *basis) {
        auto [atoms, shells] = load(natoms, anums, coords, nbasis, basis); int max_l = 0; size_t max_npgs = 0;

        for (auto shell : shells) {
            max_npgs = std::max(shell.nprim(), max_npgs); for (auto c : shell.contr) max_l = std::max(c.l, max_l);
        }

        libint2::initialize(); Engine engine(Operator::coulomb, max_npgs, max_l, 0, 1e-14);

        twoelec(I, engine, shells); libint2::finalize();
    }

    void kinetic(double *I, ulong natoms, const double *anums, const double *coords, int nbasis, const double *basis) {
        auto [atoms, shells] = load(natoms, anums, coords, nbasis, basis); int max_l = 0; size_t max_npgs = 0;

        for (auto shell : shells) {
            max_npgs = std::max(shell.nprim(), max_npgs); for (auto c : shell.contr) max_l = std::max(c.l, max_l);
        }

        libint2::initialize(); Engine engine(Operator::kinetic, max_npgs, max_l, 0, 1e-14);

        oneelec(I, engine, shells); libint2::finalize();
    }

    void nuclear(double *I, ulong natoms, const double *anums, const double *coords, int nbasis, const double *basis) {
        auto [atoms, shells] = load(natoms, anums, coords, nbasis, basis); int max_l = 0; size_t max_npgs = 0;

        for (auto shell : shells) {
            max_npgs = std::max(shell.nprim(), max_npgs); for (auto c : shell.contr) max_l = std::max(c.l, max_l);
        }

        libint2::initialize(); Engine engine(Operator::nuclear, max_npgs, max_l, 0, 1e-14);

        engine.set_params(make_point_charges(atoms));

        oneelec(I, engine, shells); libint2::finalize();
    }

    void overlap(double *I, ulong natoms, const double *anums, const double *coords, int nbasis, const double *basis) {
        auto [atoms, shells] = load(natoms, anums, coords, nbasis, basis); int max_l = 0; size_t max_npgs = 0;

        for (auto shell : shells) {
            max_npgs = std::max(shell.nprim(), max_npgs); for (auto c : shell.contr) max_l = std::max(c.l, max_l);
        }

        libint2::initialize(); Engine engine(Operator::overlap, max_npgs, max_l, 0, 1e-14);

        oneelec(I, engine, shells); libint2::finalize();
    }

    void fock(double *F, ulong natoms, const double *anums, const double *coords, ulong nbasis, const double *basis, const double *D, bool generalized) {
        auto [atoms, shells] = load(natoms, anums, coords, nbasis, basis); int max_l = 0; size_t max_npgs = 0;

        for (auto shell : shells) {
            max_npgs = std::max(shell.nprim(), max_npgs); for (auto c : shell.contr) max_l = std::max(c.l, max_l);
        }

        libint2::initialize(); Engine engine(Operator::coulomb, max_npgs, max_l, 0, 1e-14);

        generalized ? twoelec2fockg(F, engine, shells, D) : twoelec2fock(F, engine, shells, D); libint2::finalize();
    }
}
