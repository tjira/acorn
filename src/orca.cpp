#include "orca.h"

Result Orca::gradient(const System& system, const Integrals&, Result res, bool print) const {
    // define the name of the orca folder and create the gradient matrix
    std::filesystem::path orcadir = opt.folder / (".orca." + std::to_string(Timer::Now().time_since_epoch().count())); res.G = Matrix<>(system.getAtoms().size(), 3);

    // define the execution connamd and create the execution directory
    std::stringstream cmd; cmd << "./orca.sh " << system.getCharge() << " " << system.getMulti() << " " << system.getBasis(), std::filesystem::create_directory(orcadir); 

    // save the system and copy the interface
    system.save(orcadir / "molecule.xyz"), std::filesystem::copy_file(opt.interface, orcadir / "orca.sh", std::filesystem::copy_options::overwrite_existing);

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
    return res;
}
