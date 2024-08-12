#include "integral.h"
#include <argparse.hpp>
#include <filesystem>

#define FORMAT(T) [&](long ms) {char s[99]; std::sprintf(s, "%02ld:%02ld:%02ld.%03ld", ms / 3600000, ms % 3600000 / 60000, ms % 60000 / 1000, ms % 1000); return std::string(s);}(T)

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Integral Engine", "1.0", argparse::default_arguments::none);

    // define the timers
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> timers(2, std::chrono::high_resolution_clock().now()); auto tp = timers.at(0);
    auto elapsed = [](auto ms) {return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock().now() - ms).count();};

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

    // start the timer for integral calculation
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of integral calculation
    std::cout << "INTEGRALS IN AO BASIS CALCULATION: " << std::flush;

    // calculate the integrals
    libint2::initialize();
    torch::Tensor V = Integral::Nuclear(atoms, shells);
    torch::Tensor T = Integral::Kinetic(       shells);
    torch::Tensor S = Integral::Overlap(       shells);
    torch::Tensor J = Integral::Coulomb(       shells);
    libint2::finalize();

    // print the time for integral calculation
    std::cout << FORMAT(elapsed(timers.at(1))) << std::endl;

    // start the timer for integral writing
    timers.at(1) = std::chrono::high_resolution_clock().now();

    // print the header of integral writing
    std::cout << "INTEGRALS IN AO BASIS WRITING: " << std::flush;
    
    // save the integrals to disk
    torch::WriteTensor("V_AO.mat", V);
    torch::WriteTensor("T_AO.mat", T);
    torch::WriteTensor("S_AO.mat", S);
    torch::WriteTensor("J_AO.mat", J);

    // print the time for integral writing
    std::cout << FORMAT(elapsed(timers.at(1))) << std::endl;

    // print the total time
    std::cout << std::endl << "TOTAL TIME: " << FORMAT(elapsed(timers.at(0))) << std::endl;
}
