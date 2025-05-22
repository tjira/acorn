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

    void oneelec(double *ints, libint2::Engine &engine, const std::vector<Shell> &shells) {
        int nbf = 0; std::vector<size_t> sh2bf; std::vector<Engine> engines(std::getenv("OMP_NUM_THREADS") ? std::max(1, std::atoi(std::getenv("OMP_NUM_THREADS"))) : 1, engine);

        for (int i = 0, j = 0; j < shells.size(); i += shells[j].size(), j++) sh2bf.push_back(i);

        for (const auto& shell : shells) nbf += shell.size();

        #pragma omp parallel for num_threads(engines.size()) schedule(dynamic)
        for (size_t i = 0; i < shells.size(); i++) for (size_t j = i; j < shells.size(); j++) {

            int id = omp_get_thread_num(); int integral_index = 0; engines.at(id).compute(shells.at(i), shells.at(j)); if (engines.at(id).results().at(0) == nullptr) continue;

            for (size_t k = 0; k < shells.at(i).size(); k++) {
                for (size_t l = 0; l < shells.at(j).size(); l++) {
                    ints[(k + sh2bf.at(i)) * nbf + (l + sh2bf.at(j))] = engines.at(id).results().at(0)[integral_index  ];
                    ints[(l + sh2bf.at(j)) * nbf + (k + sh2bf.at(i))] = engines.at(id).results().at(0)[integral_index++];
                }
            }
        }
    }

    void twoelec(double *ints, libint2::Engine &engine, const std::vector<Shell> &shells) {
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
                            ints[(m + sh2bf.at(i)) * nbf * nbf * nbf + (n + sh2bf.at(j)) * nbf * nbf + (o + sh2bf.at(k)) * nbf + (p + sh2bf.at(l))] = engines.at(id).results().at(0)[integral_index  ];
                            ints[(m + sh2bf.at(i)) * nbf * nbf * nbf + (n + sh2bf.at(j)) * nbf * nbf + (p + sh2bf.at(l)) * nbf + (o + sh2bf.at(k))] = engines.at(id).results().at(0)[integral_index  ];
                            ints[(n + sh2bf.at(j)) * nbf * nbf * nbf + (m + sh2bf.at(i)) * nbf * nbf + (o + sh2bf.at(k)) * nbf + (p + sh2bf.at(l))] = engines.at(id).results().at(0)[integral_index  ];
                            ints[(n + sh2bf.at(j)) * nbf * nbf * nbf + (m + sh2bf.at(i)) * nbf * nbf + (p + sh2bf.at(l)) * nbf + (o + sh2bf.at(k))] = engines.at(id).results().at(0)[integral_index  ];
                            ints[(o + sh2bf.at(k)) * nbf * nbf * nbf + (p + sh2bf.at(l)) * nbf * nbf + (m + sh2bf.at(i)) * nbf + (n + sh2bf.at(j))] = engines.at(id).results().at(0)[integral_index  ];
                            ints[(o + sh2bf.at(k)) * nbf * nbf * nbf + (p + sh2bf.at(l)) * nbf * nbf + (n + sh2bf.at(j)) * nbf + (m + sh2bf.at(i))] = engines.at(id).results().at(0)[integral_index  ];
                            ints[(p + sh2bf.at(l)) * nbf * nbf * nbf + (o + sh2bf.at(k)) * nbf * nbf + (m + sh2bf.at(i)) * nbf + (n + sh2bf.at(j))] = engines.at(id).results().at(0)[integral_index  ];
                            ints[(p + sh2bf.at(l)) * nbf * nbf * nbf + (o + sh2bf.at(k)) * nbf * nbf + (n + sh2bf.at(j)) * nbf + (m + sh2bf.at(i))] = engines.at(id).results().at(0)[integral_index++];
                        }
                    }
                }
            }
        }
    }

    void coulomb(double *ints, ulong natoms, const double *anums, const double *coords, int nbasis, const double *basis) {
        auto [atoms, shells] = load(natoms, anums, coords, nbasis, basis); int max_l = 0; size_t max_npgs = 0;

        for (auto shell : shells) {
            max_npgs = std::max(shell.nprim(), max_npgs); for (auto c : shell.contr) max_l = std::max(c.l, max_l);
        }

        initialize(); Engine engine(Operator::coulomb, max_npgs, max_l, 0, 1e-14);

        twoelec(ints, engine, shells); finalize();
    }

    void kinetic(double *ints, ulong natoms, const double *anums, const double *coords, int nbasis, const double *basis) {
        auto [atoms, shells] = load(natoms, anums, coords, nbasis, basis); int max_l = 0; size_t max_npgs = 0;

        for (auto shell : shells) {
            max_npgs = std::max(shell.nprim(), max_npgs); for (auto c : shell.contr) max_l = std::max(c.l, max_l);
        }

        initialize(); Engine engine(Operator::kinetic, max_npgs, max_l, 0, 1e-14);

        oneelec(ints, engine, shells); finalize();
    }

    void nuclear(double *ints, ulong natoms, const double *anums, const double *coords, int nbasis, const double *basis) {
        auto [atoms, shells] = load(natoms, anums, coords, nbasis, basis); int max_l = 0; size_t max_npgs = 0;

        for (auto shell : shells) {
            max_npgs = std::max(shell.nprim(), max_npgs); for (auto c : shell.contr) max_l = std::max(c.l, max_l);
        }

        initialize(); Engine engine(Operator::nuclear, max_npgs, max_l, 0, 1e-14);

        engine.set_params(make_point_charges(atoms));

        oneelec(ints, engine, shells); finalize();
    }

    void overlap(double *ints, ulong natoms, const double *anums, const double *coords, int nbasis, const double *basis) {
        auto [atoms, shells] = load(natoms, anums, coords, nbasis, basis); int max_l = 0; size_t max_npgs = 0;

        for (auto shell : shells) {
            max_npgs = std::max(shell.nprim(), max_npgs); for (auto c : shell.contr) max_l = std::max(c.l, max_l);
        }

        initialize(); Engine engine(Operator::overlap, max_npgs, max_l, 0, 1e-14);

        oneelec(ints, engine, shells); finalize();
    }
}
