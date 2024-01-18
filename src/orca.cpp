#include "orca.h"

Result Orca::gradient(const System& system, const Integrals&, Result res, bool print) const {
    // define the name of the folder, where the orca calculation happens
    std::filesystem::path orcadir = opt.folder / ("orca." + std::to_string(Timer::Now().time_since_epoch().count()));

    // create the orca folder
    std::filesystem::create_directory(orcadir);

    // copy the interface to the orca folder
    std::filesystem::copy_file(opt.interface, orcadir / "orca.sh", std::filesystem::copy_options::overwrite_existing);

    // save the system to the orca folder
    system.save(orcadir / "molecule.xyz");

    // run the calculation
    auto pipe = popen(("cd " + orcadir.string() + " && ./orca.sh").c_str(), "r");

    // check for execution error
    if (!pipe) throw std::runtime_error("ORCA EXECUTION FAILED");

    // check for orca exit code
    if (pclose(pipe) != EXIT_SUCCESS) {
        throw std::runtime_error("ORCA TERMINATED UNSUCCESSFULLY");
    }

    // allocate space in the gradient matrix
    res.orca.G = Matrix<>(system.getAtoms().size(), 3);

    // extract the energy from the orca output
    std::ifstream efstream(orcadir / "energy.dat"); efstream >> res.orca.E, res.Etot = res.orca.E;

    // extract the gradient from the orca output
    std::ifstream gfstream(orcadir / "gradient.dat");
    for (size_t i = 0; i < system.getAtoms().size(); i++) {
        for (int j = 0; j < 3; j++) gfstream >> res.orca.G(i, j);
    }

    // return the results
    res.G = res.orca.G; return res;
}

Result Orca::run(const System& system, Result res, bool print) const {
    return run(system, {}, res, print);
}

Result Orca::run(const System& system, const Integrals& ints, Result res, bool print) const {
    // define the name of the folder, where the orca calculation happens
    std::filesystem::path orcadir = opt.folder / ("orca." + std::to_string(Timer::Now().time_since_epoch().count()));

    // create the orca folder
    std::filesystem::create_directory(orcadir);

    // copy the interface to the orca folder
    std::filesystem::copy_file(opt.interface, orcadir / "orca.sh", std::filesystem::copy_options::overwrite_existing);

    // save the system to the orca folder
    system.save(orcadir / "molecule.xyz");

    // run the calculation
    auto pipe = popen(("cd " + orcadir.string() + " && ./orca.sh").c_str(), "r");

    // check for execution error
    if (!pipe) throw std::runtime_error("ORCA EXECUTION FAILED");

    // check for orca exit code
    if (pclose(pipe) != EXIT_SUCCESS) {
        throw std::runtime_error("ORCA TERMINATED UNSUCCESSFULLY");
    }

    // extract the energy from the orca output
    std::ifstream fstream(orcadir / "energy.dat"); fstream >> res.orca.E, res.Etot = res.orca.E; fstream.close();

    // return the results
    return res;
}

#include "method.cpp"
template class Method<Orca>;
