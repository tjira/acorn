#include "bagel.h"

Result Bagel::gradient(const System& system, const Integrals&, Result res, bool print) const {
    // define the name of the bagel folder and create the gradient matrix
    std::filesystem::path bageldir = std::filesystem::path(ip) / (".bagel." + std::to_string(Timer::Now().time_since_epoch().count())); res.G = Matrix<>(system.getAtoms().size(), 3);

    // print the full path of the execution folder
    if (print) std::cout << "BAGEL DIRECTORY: " << bageldir << std::endl;

    // define the execution command and create the execution directory
    std::stringstream cmd; cmd << "./bagel.sh " << system.getBasis() << " " << system.getCharge() << " " << system.getMulti() << " " << opt.nstate << " " << opt.state - 1 << " > /dev/null 2>&1", std::filesystem::create_directory(bageldir);

    // save the system and copy the interface
    system.save(bageldir / "molecule.xyz"), std::filesystem::copy_file(ip / opt.interface, bageldir / "bagel.sh", std::filesystem::copy_options::overwrite_existing);

    // run the calculation and check for the error code
    auto pipe = popen(("cd " + bageldir.string() + " && " + cmd.str()).c_str(), "r");
    if (!pipe || pclose(pipe) != 0) throw std::runtime_error("EXECUTION FAILED");

    // define the energy and gradient streams
    std::ifstream gfstream(bageldir / "gradient.dat");
    std::ifstream efstream(bageldir / "energy.dat");

    // define the energies output vector
    std::vector<double> energies;

    // extract the gradient from the bagel output
    for (size_t i = 0; i < system.getAtoms().size(); i++) {
        for (int j = 0; j < 3; j++) gfstream >> res.G(i, j);
    }

    // extract the energies from the bagel output
    while (efstream >> res.Etot) {energies.push_back(res.Etot);} res.Etot = energies.at(opt.state - 1);

    // create the vector of excited energies
    res.Eexc = Vector<>(energies.size()); for (size_t i = 0; i < energies.size(); i++) res.Eexc(i) = energies.at(i);

    // return the results
    return res;
}

Result Bagel::run(const System& system, const Integrals&, Result res, bool print) const {
    // define the name of the bagel folder and create the gradient matrix
    std::filesystem::path bageldir = std::filesystem::path(ip) / (".bagel." + std::to_string(Timer::Now().time_since_epoch().count())); res.G = Matrix<>(system.getAtoms().size(), 3);

    // print the full path of the execution folder
    if (print) std::cout << "BAGEL DIRECTORY: " << bageldir << std::endl;

    // define the execution command and create the execution directory
    std::stringstream cmd; cmd << "./bagel.sh " << system.getBasis() << " " << system.getCharge() << " " << system.getMulti() << " " << opt.nstate << " " << opt.state - 1 << " > /dev/null 2>&1", std::filesystem::create_directory(bageldir);

    // save the system and copy the interface
    system.save(bageldir / "molecule.xyz"), std::filesystem::copy_file(ip / opt.interface, bageldir / "bagel.sh", std::filesystem::copy_options::overwrite_existing);

    // run the calculation and check for the error code
    auto pipe = popen(("cd " + bageldir.string() + " && " + cmd.str()).c_str(), "r");
    if (!pipe || pclose(pipe) != 0) throw std::runtime_error("EXECUTION FAILED");

    // define the energy stream and energy vector
    std::ifstream efstream(bageldir / "energy.dat"); std::vector<double> energies;;

    // extract the energies from the bagel output
    while (efstream >> res.Etot) {energies.push_back(res.Etot);} res.Etot = energies.at(opt.state - 1);

    // create the vector of excited energies
    res.Eexc = Vector<>(energies.size()); for (size_t i = 0; i < energies.size(); i++) res.Eexc(i) = energies.at(i);

    // return the results
    return res;
}
