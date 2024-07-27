#include "integral.h"
#include "timer.h"
#include <argparse.hpp>
#include <filesystem>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Integral Engine", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now(), tp;

    // add the command line arguments
    program.add_argument("-b", "--basis").help("-- Basis set used for integrals.").default_value("STO-3G");
    program.add_argument("-f", "--file").help("-- System file in the .xyz format.").default_value("molecule.xyz");
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // set the environment variable for the basis set location
    if (auto path = std::filesystem::weakly_canonical(std::filesystem::path(argv[0])).parent_path(); !std::filesystem::is_directory(std::string(DATADIR) + "/basis")) {
        #ifdef _WIN32
        _putenv_s("LIBINT_DATA_PATH", path.string().c_str());
        #else
        setenv("LIBINT_DATA_PATH", path.c_str(), true);
        #endif
    }

    // get the name of the basis and replace stars and pluses
    std::string basis = program.get("-b"); std::replace(basis.begin(), basis.end(), '*', 's'), std::replace(basis.begin(), basis.end(), '+', 'p');

    // open the system file
    std::ifstream fstream(program.get("-f")); if (!fstream.good()) throw std::runtime_error("SYSTEM FILE DOES NOT EXIST");

    // read the system file and initialize the basis set
    std::vector<libint2::Atom> atoms = libint2::read_dotxyz(fstream); libint2::BasisSet shells(basis, atoms);

    // print the timing of the system initialization
    std::cout << "SYSTEM INITIALIZATION: " << Timer::Format(Timer::Elapsed(tp)) << std::endl << std::endl;

    // calculate the integrals
    libint2::initialize();
    MEASURE("INTEGRALS IN AO BASIS CALCULATION: ",
        Matrix    V = Integral::Nuclear(atoms, shells);
        Matrix    T = Integral::Kinetic(       shells);
        Matrix    S = Integral::Overlap(       shells);
        Tensor<4> J = Integral::Coulomb(       shells);
    )
    libint2::finalize();
    
    // save the integrals to disk
    MEASURE("INTEGRALS IN AO BASIS WRITING    : ",
        Eigen::Write("V_AO.mat", V);
        Eigen::Write("T_AO.mat", T);
        Eigen::Write("S_AO.mat", S);
        Eigen::Write("J_AO.mat", J);
    )

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}
