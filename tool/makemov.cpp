// command line arguments parsing
#include <argparse/argparse.hpp>

// file operations
#include <fstream>
#include <cmath>

int main(int argc, char** argv) {
    // create the argument parser and the program timer
    argparse::ArgumentParser program("makemov", "0.2.0", argparse::default_arguments::none);

    // add the command line arguments arguments
    program.add_argument("input").help("-- Input file in .xyz format and the transformation string.").nargs(2);
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);

    // parse the command line arguments and print help if requested
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // open the input file and define the line
    std::ifstream input(program.get<std::string>("input")); std::string comment, line;

    // define the splitted transformation string and the atom groups
    std::vector<std::string> tstring; std::vector<std::vector<int>> groups;

    // check if the input file exists
    if (!input.good()) throw std::runtime_error("SYSTEM FILE DOES NOT EXIST");

    // split the transformation string
    std::istringstream iss(program.get<std::vector<std::string>>("input").at(1)); while (std::getline(iss, line, '-')) tstring.push_back(line);

    // raise an error if anything other than bond length transformation is requested
    if (tstring.at(0) != "b") throw std::runtime_error("UNKNOWN TRANSFORMATION REQUESTED");

    // extract the atom groups from the transformation string
    for (size_t i = 1; i < tstring.size() - 3; i++) {
        std::vector<int> group; std::istringstream iss(tstring.at(i)); while (std::getline(iss, line, ',')) {group.push_back(std::stoi(line));}; groups.push_back(group);
    }

    // extract the range and number of steps from the transformation string
    double start = std::stod(tstring.at(tstring.size() - 3)), end = std::stod(tstring.at(tstring.size() - 2)); int steps = std::stod(tstring.at(tstring.size() - 1));

    while (std::getline(input, line)) {
        // read the comment and convert the number of atoms to integer
        std::getline(input, comment); int natoms = std::stoi(line);

        // define the atomic symbols and coordinates
        std::vector<std::string> symbols(natoms); std::vector<double> x(natoms), y(natoms), z(natoms);

        // read the atomic symbols and coordinates
        for (int i = 0; i < natoms; i++) {
            std::getline(input, line); std::istringstream iss(line); iss >> symbols.at(i) >> x.at(i) >> y.at(i) >> z.at(i);
        }

        // bond length transformation
        if (tstring.at(0) == "b") {
            // define the coordinates of the first atom of the bond
            double x1 = x.at(groups.at(0).at(groups.at(0).size() - 1) - 1), y1 = y.at(groups.at(0).at(groups.at(0).size() - 1) - 1), z1 = z.at(groups.at(0).at(groups.at(0).size() - 1) - 1);

            // define the coordinates of the second atom of the bond
            double x2 = x.at(groups.at(1).at(0) - 1), y2 = y.at(groups.at(1).at(0) - 1), z2 = z.at(groups.at(1).at(0) - 1);

            // get the vector of the bond and its length
            double dx = x2 - x1, dy = y2 - y1, dz = z2 - z1, n = std::sqrt(dx * dx + dy * dy + dz * dz);

            // perform the bond length transformation
            for (int i = 0; i < steps; i++) {
                // define the additive factors for this iteration
                double fx = dx / n * (start + i * (end - start) / (steps - 1)) - dx;
                double fy = dy / n * (start + i * (end - start) / (steps - 1)) - dy;
                double fz = dz / n * (start + i * (end - start) / (steps - 1)) - dz;

                // print the number of atoms and comment
                std::printf("%d\n%s\n", natoms, comment.c_str());

                // print the centered coordinates
                for (int j = 0; j < natoms; j++) {
                    if(std::find(groups.at(0).begin(), groups.at(0).end(), j + 1) == groups.at(0).end()) {
                        std::printf("%2s %20.14f %20.14f %20.14f\n", symbols.at(j).c_str(), x.at(j) + fx, y.at(j) + fy, z.at(j) + fz);
                    } else std::printf("%2s %20.14f %20.14f %20.14f\n", symbols.at(j).c_str(), x.at(j), y.at(j), z.at(j));
                }
            }
        }
    }
}
