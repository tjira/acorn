#define DATADIR "/home/tojira/Projects/acorn"

#include <libint2.hpp>

#include <fstream>

extern "C" {
    using namespace libint2;

    void oneelec(double *ints, libint2::Engine &engine, const BasisSet &obs) {
        int nbf = obs.nbf(); std::vector<size_t> sh2bf = obs.shell2bf();

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

    void twoelec(double *ints, libint2::Engine &engine, const BasisSet &obs) {
        int nbf = obs.nbf(); std::vector<size_t> sh2bf = obs.shell2bf();

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

    void coulomb(double *ints, const char *system, const char *basis) {
        std::ifstream system_stream(system); std::vector<Atom> atoms = read_dotxyz(system_stream); BasisSet obs(basis, atoms);

        initialize(); Engine engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 0, 1e-14);

        twoelec(ints, engine, obs); finalize();
    }

    void kinetic(double *ints, const char *system, const char *basis) {
        std::ifstream system_stream(system); std::vector<Atom> atoms = read_dotxyz(system_stream); BasisSet obs(basis, atoms);

        initialize(); Engine engine(Operator::kinetic, obs.max_nprim(), obs.max_l(), 0, 1e-14);

        oneelec(ints, engine, obs); finalize();
    }

    void nuclear(double *ints, const char *system, const char *basis) {
        std::ifstream system_stream(system); std::vector<Atom> atoms = read_dotxyz(system_stream); BasisSet obs(basis, atoms);

        initialize(); Engine engine(Operator::nuclear, obs.max_nprim(), obs.max_l(), 0, 1e-14);

        engine.set_params(make_point_charges(atoms));

        oneelec(ints, engine, obs); finalize();
    }

    void overlap(double *ints, const char *system, const char *basis) {
        std::ifstream system_stream(system); std::vector<Atom> atoms = read_dotxyz(system_stream); BasisSet obs(basis, atoms);

        initialize(); Engine engine(Operator::overlap, obs.max_nprim(), obs.max_l(), 0, 1e-14);

        oneelec(ints, engine, obs); finalize();
    }
}
