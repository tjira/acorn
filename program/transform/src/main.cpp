#include "transform.h"
#include "timer.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Hartree-Fock Integral Transform Engine", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now(), tp;

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
        MEASURE("INTEGRALS IN AO BASIS & COEFFICIENT MATRIX IN MO BASIS READING: ",
            E = Eigen::LoadMatrix("E_MO.mat");
            C = Eigen::LoadMatrix("C_MO.mat");
            V = Eigen::LoadMatrix("V_AO.mat");
            T = Eigen::LoadMatrix("T_AO.mat");
            S = Eigen::LoadMatrix("S_AO.mat");
            J = Eigen::LoadTensor("J_AO.mat");
            F = Eigen::LoadMatrix("F_AO.mat");
        )
    }

    // print a new line
    if (program.get<bool>("-o") || program.get<bool>("-s")) std::cout << std::endl;

    // integrals in MO basis calculation
    if (program.get<bool>("-o")) {
        MEASURE("INTEGRAL TRANSFORMS TO MO BASIS: ",
            Jmo = Transform::CoulombSpatial(J, C);
            Vmo = Transform::SingleSpatial (V, C);
            Tmo = Transform::SingleSpatial (T, C);
            Smo = Transform::SingleSpatial (S, C);
            Fmo = Transform::SingleSpatial (F, C);
        )
    }

    // integrals in MO basis calculation
    if (program.get<bool>("-s")) {
        MEASURE("INTEGRAL TRANSFORMS TO MS BASIS: ",
            Jms = Transform::CoulombSpin(J, C);
            Vms = Transform::SingleSpin (V, C);
            Tms = Transform::SingleSpin (T, C);
            Sms = Transform::SingleSpin (S, C);
            Fms = Transform::SingleSpin (F, C);
            Ems = Vector::NullaryExpr(2 * E.rows(), [&](int i){return E(i / 2);});
        )
    }

    // integrals in molecular spatial orbital basis writing
    if (program.get<bool>("-o")) {
        MEASURE("INTEGRALS IN MO BASIS WRITING:   ",
            Eigen::Write("V_MO.mat", Vmo);
            Eigen::Write("T_MO.mat", Tmo);
            Eigen::Write("S_MO.mat", Smo);
            Eigen::Write("J_MO.mat", Jmo);
            Eigen::Write("F_MO.mat", Fmo);
        )
    }

    // integrals in molecular spatial orbital basis writing
    if (program.get<bool>("-s")) {
        MEASURE("INTEGRALS IN MS BASIS WRITING:   ",
            Eigen::Write("V_MS.mat", Vms);
            Eigen::Write("T_MS.mat", Tms);
            Eigen::Write("S_MS.mat", Sms);
            Eigen::Write("J_MS.mat", Jms);
            Eigen::Write("E_MS.mat", Ems);
            Eigen::Write("F_MS.mat", Fms);
        )
    }

    // print the total time
    std::cout << (program.get<bool>("-o") || program.get<bool>("-s") ? "\n" : "") << "TOTAL TIME: " << Timer::Format(Timer::Elapsed(start)) << std::endl;
}
