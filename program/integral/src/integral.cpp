#include "integral.h"

torch::Tensor Acorn::Integral::Single(libint2::Engine& engine, const libint2::BasisSet& shells) {
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

torch::Tensor Acorn::Integral::Double(libint2::Engine& engine, const libint2::BasisSet& shells) {
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

torch::Tensor Acorn::Integral::Nuclear(const std::vector<libint2::Atom>& atoms, const libint2::BasisSet& shells) {
    libint2::Engine engine(libint2::Operator::nuclear, shells.max_nprim(), shells.max_l(), 0, 1e-12);
    engine.set_params(libint2::make_point_charges(atoms)); return Single(engine, shells);
}

torch::Tensor Acorn::Integral::Kinetic(const libint2::BasisSet& shells) {
    libint2::Engine engine(libint2::Operator::kinetic, shells.max_nprim(), shells.max_l(), 0, 1e-12);
    return Single(engine, shells);
}

torch::Tensor Acorn::Integral::Overlap(const libint2::BasisSet& shells) {
    libint2::Engine engine(libint2::Operator::overlap, shells.max_nprim(), shells.max_l(), 0, 1e-12);
    return Single(engine, shells);
}

torch::Tensor Acorn::Integral::Coulomb(const libint2::BasisSet& shells) {
    libint2::Engine engine(libint2::Operator::coulomb, shells.max_nprim(), shells.max_l(), 0, 1e-12);
    return Double(engine, shells);
}

void Acorn::Integral::run(const Options& opt, std::vector<timepoint>& timers) {
    // set the environment variable for the basis set location
    if (auto path = std::filesystem::weakly_canonical(std::filesystem::path(opt.executable)).parent_path(); !std::filesystem::is_directory(std::string(DATADIR) + "/basis")) {
        #ifdef _WIN32
        _putenv_s("LIBINT_DATA_PATH", path.string().c_str());
        #else
        setenv("LIBINT_DATA_PATH", path.c_str(), true);
        #endif
    }

    // open the system file
    std::ifstream fstream(opt.file); if (!fstream.good()) throw std::runtime_error("SYSTEM FILE DOES NOT EXIST");

    // read the system file and initialize the basis set
    std::vector<libint2::Atom> atoms = libint2::read_dotxyz(fstream); libint2::BasisSet shells(opt.basis, atoms);

    // start the timer for integral calculation
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of integral calculation
    std::cout << "INTEGRALS IN AO BASIS CALCULATION: " << std::flush;

    // calculate the integrals
    libint2::initialize();
    torch::Tensor V = Acorn::Integral::Nuclear(atoms, shells);
    torch::Tensor T = Acorn::Integral::Kinetic(       shells);
    torch::Tensor S = Acorn::Integral::Overlap(       shells);
    torch::Tensor J = Acorn::Integral::Coulomb(       shells);
    libint2::finalize();

    // print the time for integral calculation
    std::cout << eltime(timers.at(1)) << std::endl;

    // start the timer for integral writing
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of integral writing
    std::cout << "INTEGRALS IN AO BASIS WRITING:     " << std::flush;
    
    // save the integrals to disk
    torch::WriteTensor("V_AO.mat", V);
    torch::WriteTensor("T_AO.mat", T);
    torch::WriteTensor("S_AO.mat", S);
    torch::WriteTensor("J_AO.mat", J);

    // print the time for integral writing
    std::cout << eltime(timers.at(1)) << std::endl;

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << eltime(timers.at(0)) << std::endl;
}
