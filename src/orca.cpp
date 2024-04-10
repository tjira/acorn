#include "orca.h"

Result Orca::gradient(const System& system, const Integrals&, Result res, bool print) const {
    // define the name of the orca folder and create the gradient matrix
    std::filesystem::path orcadir = std::filesystem::path(ip) / (".orca." + std::to_string(Timer::Now().time_since_epoch().count())); res.G = Matrix<>(system.getAtoms().size(), 3);

    // print the full path of the execution folder
    if (print) std::cout << "ORCA DIRECTORY: " << orcadir << std::endl;

    // define the execution command
    std::stringstream cmd; cmd << "./orca.sh " << system.getBasis() << " " << system.getCharge() << " " << system.getMulti() << " " << opt.method << " > /dev/null 2>&1";

    // create the execution directory
    std::filesystem::create_directory(orcadir);

    // save the system and copy the interface
    system.save((orcadir / "molecule.xyz").string()), std::filesystem::copy_file(ip / opt.interface, orcadir / "orca.sh", std::filesystem::copy_options::overwrite_existing);

    // run the calculation and check for the error code
    auto pipe = popen(("cd " + orcadir.string() + " && " + cmd.str()).c_str(), "r");
    if (!pipe || pclose(pipe) != 0) throw std::runtime_error("EXECUTION FAILED");

    // define the energy and gradient streams
    std::ifstream gfstream(orcadir / "gradient.dat");
    std::ifstream efstream(orcadir / "energy.dat");

    // extract the gradient and energy from the orca output
    for (size_t i = 0; i < system.getAtoms().size(); i++) {
        for (int j = 0; j < 3; j++) gfstream >> res.G(i, j);
    } efstream >> res.Etot;

    // return the results
    res.Eexc = Vector<>(1); res.Eexc << res.Etot; return res;
}

Result Orca::run(const System& system, const Integrals&, Result res, bool print) const {
    // define the name of the orca folder and create the gradient matrix
    std::filesystem::path orcadir = std::filesystem::path(ip) / (".orca." + std::to_string(Timer::Now().time_since_epoch().count())); res.G = Matrix<>(system.getAtoms().size(), 3);

    // print the full path of the execution folder
    if (print) std::cout << "ORCA DIRECTORY: " << orcadir << std::endl;

    // define the execution command
    std::stringstream cmd; cmd << "./orca.sh " << system.getBasis() << " " << system.getCharge() << " " << system.getMulti() << " " << opt.method << " > /dev/null 2>&1";

    // create the execution directory
    std::filesystem::create_directory(orcadir);

    // save the system and copy the interface
    system.save((orcadir / "molecule.xyz").string()), std::filesystem::copy_file(ip / opt.interface, orcadir / "orca.sh", std::filesystem::copy_options::overwrite_existing);

    // run the calculation and check for the error code
    auto pipe = popen(("cd " + orcadir.string() + " && " + cmd.str()).c_str(), "r");
    if (!pipe || pclose(pipe) != 0) throw std::runtime_error("EXECUTION FAILED");

    // define the energy and gradient streams
    std::ifstream efstream(orcadir / "energy.dat"); efstream >> res.Etot;

    // return the results
    res.Eexc = Vector<>(1); res.Eexc << res.Etot; return res;
}
