#define DATADIR "/home/tojira/Projects/acorn"

#include <libint2.hpp>

#include <fstream>

size_t max_nprim(const std::vector<libint2::Shell>& shells) {
  size_t n = 0;
  for (auto shell : shells) n = std::max(shell.nprim(), n);
  return n;
}

int max_l(const std::vector<libint2::Shell>& shells) {
  int l = 0;
  for (auto shell : shells)
    for (auto c : shell.contr) l = std::max(c.l, l);
  return l;
}

std::vector<size_t> shell2bf(const std::vector<libint2::Shell>& shells) {
    std::vector<size_t> result; result.reserve(shells.size());

    size_t n = 0;
    for (auto shell : shells) {
        result.push_back(n);
        n += shell.size();
    }

    return result;
}

extern "C" {
    using namespace libint2; typedef double real_t;

    void oneelec(double *ints, libint2::Engine &engine, const std::vector<Shell> &obs) {
        int nbf = 0; std::vector<size_t> sh2bf = shell2bf(obs);

        for (const auto& shell : obs) nbf += shell.size();

        for (size_t i = 0; i < obs.size(); i++) for (size_t j = i; j < obs.size(); j++) {

            int integral_index = 0; engine.compute(obs.at(i), obs.at(j)); if (engine.results().at(0) == nullptr) continue;

            for (size_t k = 0; k < obs.at(i).size(); k++) {
                for (size_t l = 0; l < obs.at(j).size(); l++) {
                    ints[(k + sh2bf.at(i)) * nbf + (l + sh2bf.at(j))] = engine.results().at(0)[integral_index  ];
                    ints[(l + sh2bf.at(j)) * nbf + (k + sh2bf.at(i))] = engine.results().at(0)[integral_index++];
                }
            }
        }
    }

    void twoelec(double *ints, libint2::Engine &engine, const std::vector<Shell> &obs) {
        int nbf = 0; std::vector<size_t> sh2bf = shell2bf(obs);

        for (const auto& shell : obs) nbf += shell.size();

        for (size_t i = 0; i < obs.size(); i++) for (size_t j = i; j < obs.size(); j++) for (size_t k = i; k < obs.size(); k++) for (size_t l = (i == k ? j : k); l < obs.size(); l++) {

            int integral_index = 0; engine.compute(obs.at(i), obs.at(j), obs.at(k), obs.at(l)); if (engine.results().at(0) == nullptr) continue;

            for (size_t m = 0; m < obs.at(i).size(); m++) {
                for (size_t n = 0; n < obs.at(j).size(); n++) {
                    for (size_t o = 0; o < obs.at(k).size(); o++) {
                        for (size_t p = 0; p < obs.at(l).size(); p++) {
                            ints[(m + sh2bf.at(i)) * nbf * nbf * nbf + (n + sh2bf.at(j)) * nbf * nbf + (o + sh2bf.at(k)) * nbf + (p + sh2bf.at(l))] = engine.results().at(0)[integral_index  ];
                            ints[(m + sh2bf.at(i)) * nbf * nbf * nbf + (n + sh2bf.at(j)) * nbf * nbf + (p + sh2bf.at(l)) * nbf + (o + sh2bf.at(k))] = engine.results().at(0)[integral_index  ];
                            ints[(n + sh2bf.at(j)) * nbf * nbf * nbf + (m + sh2bf.at(i)) * nbf * nbf + (o + sh2bf.at(k)) * nbf + (p + sh2bf.at(l))] = engine.results().at(0)[integral_index  ];
                            ints[(n + sh2bf.at(j)) * nbf * nbf * nbf + (m + sh2bf.at(i)) * nbf * nbf + (p + sh2bf.at(l)) * nbf + (o + sh2bf.at(k))] = engine.results().at(0)[integral_index  ];
                            ints[(o + sh2bf.at(k)) * nbf * nbf * nbf + (p + sh2bf.at(l)) * nbf * nbf + (m + sh2bf.at(i)) * nbf + (n + sh2bf.at(j))] = engine.results().at(0)[integral_index  ];
                            ints[(o + sh2bf.at(k)) * nbf * nbf * nbf + (p + sh2bf.at(l)) * nbf * nbf + (n + sh2bf.at(j)) * nbf + (m + sh2bf.at(i))] = engine.results().at(0)[integral_index  ];
                            ints[(p + sh2bf.at(l)) * nbf * nbf * nbf + (o + sh2bf.at(k)) * nbf * nbf + (m + sh2bf.at(i)) * nbf + (n + sh2bf.at(j))] = engine.results().at(0)[integral_index  ];
                            ints[(p + sh2bf.at(l)) * nbf * nbf * nbf + (o + sh2bf.at(k)) * nbf * nbf + (n + sh2bf.at(j)) * nbf + (m + sh2bf.at(i))] = engine.results().at(0)[integral_index++];
                        }
                    }
                }
            }
        }
    }

    void coulomb(double *ints, int natoms, const double *anums, const double *coords, int nbasis, const double *basis) {
        std::vector<Atom> atoms(natoms); std::vector<Shell> obs;

        for (int i = 0; i < natoms; i++) {
            atoms.at(i) = {(int)anums[i], coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]};
        }

        for (int i = 0; i < nbasis;) {

            int n = basis[i]; int am = basis[i + 1]; double x = basis[i + 2]; double y = basis[i + 3]; double z = basis[i + 4];

            svector<double> alpha(n), c(n);

            for (int j = 0; j < n; j++) {alpha[j] = basis[i + 5 + j]; c[j] = basis[i + 5 + n + j];}

            obs.push_back(Shell{alpha, {{am, false, c}}, {{x, y, z}}}); i += 2 * n + 5;
        }

        initialize(); Engine engine(Operator::coulomb, max_nprim(obs), max_l(obs), 0, 1e-14);

        twoelec(ints, engine, obs); finalize();
    }

    void kinetic(double *ints, int natoms, const double *anums, const double *coords, int nbasis, const double *basis) {
        std::vector<Atom> atoms(natoms); std::vector<Shell> obs;

        for (int i = 0; i < natoms; i++) {
            atoms.at(i) = {(int)anums[i], coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]};
        }

        for (int i = 0; i < nbasis;) {

            int n = basis[i]; int am = basis[i + 1]; double x = basis[i + 2]; double y = basis[i + 3]; double z = basis[i + 4];

            svector<double> alpha(n), c(n);

            for (int j = 0; j < n; j++) {alpha[j] = basis[i + 5 + j]; c[j] = basis[i + 5 + n + j];}

            obs.push_back(Shell{alpha, {{am, false, c}}, {{x, y, z}}}); i += 2 * n + 5;
        }

        initialize(); Engine engine(Operator::kinetic, max_nprim(obs), max_l(obs), 0, 1e-14);

        oneelec(ints, engine, obs); finalize();
    }

    void nuclear(double *ints, int natoms, const double *anums, const double *coords, int nbasis, const double *basis) {
        std::vector<Atom> atoms(natoms); std::vector<Shell> obs;

        for (int i = 0; i < natoms; i++) {
            atoms.at(i) = {(int)anums[i], coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]};
        }

        for (int i = 0; i < nbasis;) {

            int n = basis[i]; int am = basis[i + 1]; double x = basis[i + 2]; double y = basis[i + 3]; double z = basis[i + 4];

            svector<double> alpha(n), c(n);

            for (int j = 0; j < n; j++) {alpha[j] = basis[i + 5 + j]; c[j] = basis[i + 5 + n + j];}

            obs.push_back(Shell{alpha, {{am, false, c}}, {{x, y, z}}}); i += 2 * n + 5;
        }

        initialize(); Engine engine(Operator::nuclear, max_nprim(obs), max_l(obs), 0, 1e-14);

        engine.set_params(make_point_charges(atoms));

        oneelec(ints, engine, obs); finalize();
    }

    void overlap(double *ints, int natoms, const double *anums, const double *coords, int nbasis, const double *basis) {
        std::vector<Atom> atoms(natoms); std::vector<Shell> obs;

        for (int i = 0; i < natoms; i++) {
            atoms.at(i) = {(int)anums[i], coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]};
        }

        for (int i = 0; i < nbasis;) {

            int n = basis[i]; int am = basis[i + 1]; double x = basis[i + 2]; double y = basis[i + 3]; double z = basis[i + 4];

            svector<double> alpha(n), c(n);

            for (int j = 0; j < n; j++) {
                alpha[j] = basis[i + 5 + j]; c[j] = basis[i + 5 + n + j];
            }

            obs.push_back(Shell{alpha, {{am, false, c}}, {{x, y, z}}}); i += 2 * n + 5;
        }

        initialize(); Engine engine(Operator::overlap, max_nprim(obs), max_l(obs), 0, 1e-14);

        oneelec(ints, engine, obs); finalize();
    }
}
