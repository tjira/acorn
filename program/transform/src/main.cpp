#include "transform.h"
#include <argparse.hpp>

#define FORMAT(T) [&](long ms) {char s[99]; std::sprintf(s, "%02ld:%02ld:%02ld.%03ld", ms / 3600000, ms % 3600000 / 60000, ms % 60000 / 1000, ms % 1000); return std::string(s);}(T)

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Integral Transform Engine", "1.0", argparse::default_arguments::none);

    // define the timers
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> timers(2, std::chrono::high_resolution_clock().now()); auto tp = timers.at(0);
    auto elapsed = [](auto ms) {return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock().now() - ms).count();};

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-o", "--orbital").help("-- Transforms to the basis of spatial molecular orbitals.").default_value(false).implicit_value(true);
    program.add_argument("-s", "--spinorbital").help("-- Transforms to the basis of molecular spinorbitals.").default_value(false).implicit_value(true);

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}

    // define all the matrices and tensors used throughout the program
    Matrix C, E, V, T, S, F, Vmo, Tmo, Smo, Fmo, Ems, Vms, Tms, Sms, Fms; Tensor<4> J, Jmo, Jms;

    // load the coefficient matrix and integrals in AO basis from disk
    if (program.get<bool>("-o") || program.get<bool>("-s")) {

        // start the timer for integral loading
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // print the header of integral loading
        std::cout << "SYSTEM AND INTEGRALS IN AO BASIS READING: " << std::flush;

        // load the integrals in AO basis
        E = Eigen::LoadVector("E_MO.mat");
        C = Eigen::LoadMatrix("C_MO.mat");
        V = Eigen::LoadMatrix("V_AO.mat");
        T = Eigen::LoadMatrix("T_AO.mat");
        S = Eigen::LoadMatrix("S_AO.mat");
        J = Eigen::LoadTensor("J_AO.mat");
        F = Eigen::LoadMatrix("F_AO.mat");

        // print the time for integral loading
        std::cout << FORMAT(elapsed(timers.at(1))) << std::endl;
    }

    // print a new line
    if (program.get<bool>("-o") || program.get<bool>("-s")) std::cout << std::endl;

    // integrals in MO basis calculation
    if (program.get<bool>("-o")) {
        
        // start the timer for integral transformation
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // print the header of integral transformation
        std::cout << "INTEGRAL TRANSFORMS TO MO BASIS: " << std::flush;

        // calculate the integrals in MO basis
        Jmo = Transform::CoulombSpatial(J, C);
        Vmo = Transform::SingleSpatial (V, C);
        Tmo = Transform::SingleSpatial (T, C);
        Smo = Transform::SingleSpatial (S, C);
        Fmo = Transform::SingleSpatial (F, C);

        // print the time for integral transformation
        std::cout << FORMAT(elapsed(timers.at(1))) << std::endl;
    }

    // integrals in MO basis calculation
    if (program.get<bool>("-s")) {
    
        // start the timer for integral transformation
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // print the header of integral transformation
        std::cout << "INTEGRAL TRANSFORMS TO MS BASIS: " << std::flush;

        // calculate the integrals in MS basis
        Jms = Transform::CoulombSpin(J, C);
        Vms = Transform::SingleSpin (V, C);
        Tms = Transform::SingleSpin (T, C);
        Sms = Transform::SingleSpin (S, C);
        Fms = Transform::SingleSpin (F, C);
        Ems = Vector::NullaryExpr(2 * E.rows(), [&](int i){return E(i / 2);});

        // print the time for integral transformation
        std::cout << FORMAT(elapsed(timers.at(1))) << std::endl;
    }

    // integrals in molecular spatial orbital basis writing
    if (program.get<bool>("-o")) {

        // start the timer for integral writing
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // print the header of integral writing
        std::cout << "INTEGRALS IN MO BASIS WRITING:   " << std::flush;

        // save the integrals to disk
        Eigen::Write("V_MO.mat", Vmo);
        Eigen::Write("T_MO.mat", Tmo);
        Eigen::Write("S_MO.mat", Smo);
        Eigen::Write("J_MO.mat", Jmo);
        Eigen::Write("F_MO.mat", Fmo);

        // print the time for integral writing
    }

    // integrals in molecular spatial orbital basis writing
    if (program.get<bool>("-s")) {

        // start the timer for integral writing
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // print the header of integral writing
        std::cout << "INTEGRALS IN MS BASIS WRITING:   " << std::flush;

        // save the integrals to disk
        Eigen::Write("V_MS.mat", Vms);
        Eigen::Write("T_MS.mat", Tms);
        Eigen::Write("S_MS.mat", Sms);
        Eigen::Write("J_MS.mat", Jms);
        Eigen::Write("E_MS.mat", Ems);
        Eigen::Write("F_MS.mat", Fms);

        // print the time for integral writing
        std::cout << FORMAT(elapsed(timers.at(1))) << std::endl;
    }

    // print the total time
    std::cout << (program.get<bool>("-o") || program.get<bool>("-s") ? "\n" : "") << "TOTAL TIME: " << FORMAT(elapsed(timers.at(0))) << std::endl;
}
