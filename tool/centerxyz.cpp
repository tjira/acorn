// command line arguments parsing
#include <argparse/argparse.hpp>

// file operations
#include <fstream>

int main(int argc, char** argv) {
    // create the argument parser and the program timer
    argparse::ArgumentParser program("centerxyz", "0.2.0", argparse::default_arguments::none);

    // add the command line arguments arguments
    program.add_argument("input").help("-- Input file in .xyz format to center.");
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);

    // parse the command line arguments and print help if requested
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // open the input file and define line
    std::ifstream input(program.get<std::string>("input")); std::string comment, line;

    // check if the input file exists
    if (!input.good()) throw std::runtime_error("SYSTEM FILE DOES NOT EXIST");

    while (std::getline(input, line)) {
        // read the comment and convert the number of atoms to integer
        std::getline(input, comment); int natoms = std::stoi(line);

        // define the atomic symbols and coordinates
        std::vector<std::string> symbols(natoms); std::vector<double> x(natoms), y(natoms), z(natoms);

        // read the atomic symbols and coordinates
        for (int i = 0; i < natoms; i++) {
            std::getline(input, line); std::istringstream iss(line); iss >> symbols.at(i) >> x.at(i) >> y.at(i) >> z.at(i);
        }

        // define the geometric center of the molecule
        double xcenter = 0.0, ycenter = 0.0, zcenter = 0.0;

        // calculate the geometric center of the molecule
        for (int i = 0; i < natoms; i++) {
            xcenter += x.at(i); ycenter += y.at(i); zcenter += z.at(i);
        }

        // divide by the number of atoms to get the center
        xcenter /= (double)natoms; ycenter /= (double)natoms; zcenter /= (double)natoms;

        // print the number of atoms and comment
        std::printf("%d\n%s\n", natoms, comment.c_str());

        // print the centered coordinates
        for (int i = 0; i < natoms; i++) {
            std::printf("%2s %20.14f %20.14f %20.14f\n", symbols.at(i).c_str(), x.at(i) - xcenter, y.at(i) - ycenter, z.at(i) - zcenter);
        }
    }
}
