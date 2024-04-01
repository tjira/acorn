#include "bagel.h"

Result Bagel::gradient(const System& system, const Integrals&, Result res, bool print) const {
    // define the name of the bagel folder and create the gradient matrix
    std::filesystem::path bageldir = std::filesystem::path(ip) / (".bagel." + std::to_string(Timer::Now().time_since_epoch().count())); res.G = Matrix<>(system.getAtoms().size(), 3);

    // print the full path of the execution folder
    if (print) std::cout << "BAGEL DIRECTORY: " << bageldir << std::endl;

    // define the execution command and create the execution directory
    std::stringstream cmd; cmd << "./bagel.sh " << system.getBasis() << " " << system.getCharge() << " " << system.getMulti() << " > /dev/null 2>&1", std::filesystem::create_directory(bageldir);

    // save the system and copy the interface
    system.save(bageldir / "molecule.xyz"), std::filesystem::copy_file(ip / opt.interface, bageldir / "bagel.sh", std::filesystem::copy_options::overwrite_existing);

    // run the calculation and check for the error code
    auto pipe = popen(("cd " + bageldir.string() + " && " + cmd.str()).c_str(), "r");
    if (!pipe || pclose(pipe) != 0) throw std::runtime_error("EXECUTION FAILED");

    // define the energy and gradient streams
    std::ifstream gfstream(bageldir / "gradient.dat");
    std::ifstream efstream(bageldir / "energy.dat");

    // extract the gradient and energy from the bagel output
    for (size_t i = 0; i < system.getAtoms().size(); i++) {
        for (int j = 0; j < 3; j++) gfstream >> res.G(i, j);
    } efstream >> res.Etot;

    // return the results
    res.Eexc = Vector<>(1); res.Eexc << res.Etot; return res;
}

Result Bagel::run(const System& system, const Integrals&, Result res, bool print) const {
    // define the name of the bagel folder and create the gradient matrix
    std::filesystem::path bageldir = std::filesystem::path(ip) / (".bagel." + std::to_string(Timer::Now().time_since_epoch().count())); res.G = Matrix<>(system.getAtoms().size(), 3);

    // print the full path of the execution folder
    if (print) std::cout << "BAGEL DIRECTORY: " << bageldir << std::endl;

    // define the execution command and create the execution directory
    std::stringstream cmd; cmd << "./bagel.sh " << system.getBasis() << " " << system.getCharge() << " " << system.getMulti() << " > /dev/null 2>&1", std::filesystem::create_directory(bageldir);

    // save the system and copy the interface
    system.save(bageldir / "molecule.xyz"), std::filesystem::copy_file(ip / opt.interface, bageldir / "bagel.sh", std::filesystem::copy_options::overwrite_existing);

    // run the calculation and check for the error code
    auto pipe = popen(("cd " + bageldir.string() + " && " + cmd.str()).c_str(), "r");
    if (!pipe || pclose(pipe) != 0) throw std::runtime_error("EXECUTION FAILED");

    // define the energy and gradient streams
    std::ifstream efstream(bageldir / "energy.dat"); efstream >> res.Etot;

    // return the results
    res.Eexc = Vector<>(1); res.Eexc << res.Etot; return res;
}
