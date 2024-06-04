#include "system.h"
#include "timer.h"
#include <argparse.hpp>

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Acorn Moller-Plesset Pertrubation Theory Program", "1.0", argparse::default_arguments::none); Timer::Timepoint start = Timer::Now();

    // add the command line arguments
    program.add_argument("-h", "--help").help("-- This help message.").default_value(false).implicit_value(true);
    program.add_argument("-f", "--file").help("-- System file in the .xyz format.").default_value("molecule.xyz");
    program.add_argument("-o", "--order").help("-- Order of theory.").default_value(2).scan<'i', int>();

    // parse the command line arguments
    try {program.parse_args(argc, argv);} catch (const std::runtime_error& error) {
        if (!program.get<bool>("-h")) {std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);}
    } if (program.get<bool>("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);} Timer::Timepoint tp = Timer::Now();

    // load the system and integrals in MS basis from disk
    MEASURE("SYSTEM AND INTEGRALS IN MS BASIS READING: ",
        Matrix    Vms = Eigen::LoadMatrix("V_MS.mat");
        Matrix    Tms = Eigen::LoadMatrix("T_MS.mat");
        Tensor<4> Jms = Eigen::LoadTensor("J_MS.mat");
        Matrix    Ems = Eigen::LoadMatrix("E_MS.mat");
        System system(program.get("-f"));
    )

    // extract the number of occupied and virtual orbitals and define the energy
    int nocc = system.nocc(); int nvirt = Jms.dimension(0) / 2 - nocc; double E = 0; 

    // initialize the antisymmetrized Coulomb integrals and the Hamiltonian matrix in MS basis
    Tensor<4> Jmsa = Jms - Jms.shuffle(Eigen::array<int, 4>{0, 3, 2, 1}); Matrix Hms = Tms + Vms;

     // define the Coulomb integrals slice indices
    Eigen::array<Eigen::Index, 4> vvvvis = {2 * nocc, 2 * nocc, 2 * nocc, 2 * nocc}, vvvvie = {2 * nvirt, 2 * nvirt, 2 * nvirt, 2 * nvirt};
    Eigen::array<Eigen::Index, 4> ovovis = {0, 2 * nocc, 0, 2 * nocc}, ovovie = {2 * nocc, 2 * nvirt, 2 * nocc, 2 * nvirt};
    Eigen::array<Eigen::Index, 4> oovvis = {0, 0, 2 * nocc, 2 * nocc}, oovvie = {2 * nocc, 2 * nocc, 2 * nvirt, 2 * nvirt};
    Eigen::array<Eigen::Index, 4> oooois = {0, 0, 0, 0}, ooooie = {2 * nocc, 2 * nocc, 2 * nocc, 2 * nocc};

    // obtain the Coulomb integral slices and define the amplitude tensor
    Tensor<4> oooo = Jmsa.slice(oooois, ooooie), oovv = Jmsa.slice(oovvis, oovvie), ovov = Jmsa.slice(ovovis, ovovie), vvvv = Jmsa.slice(vvvvis, vvvvie); Tensor<4> t = ovov;

    // divide the amplitudes by the orbital energy differences
    for (int i = 0; i < 2 * nocc; i++) {
        for (int j = 0; j < 2 * nocc; j++) {
            for (int a = 0; a < 2 * nvirt; a++) {
                for (int b = 0; b < 2 * nvirt; b++) {
                    t(i, a, j, b) /= Ems(2 * nocc + a) + Ems(2 * nocc + b) - Ems(i) - Ems(j);
                }
            }
        }
    }

    // calculate the HF energy
    for (int i = 0; i < 2 * nocc; i++) {
        E += Hms(i, i); for (int j = 0; j < 2 * nocc; j++) E += 0.5 * Jmsa(i, i, j, j);
    }

    // print new line
    std::cout << std::endl;

    // add the MP2 correlation energy
    if (tp = Timer::Now(); program.get<int>("-o") > 1) {
        // define the contraction axes
        Eigen::IndexPair<int> a1(0, 0), a2(1, 1), a3(2, 2), a4(3, 3);

        // add the MP2 correlation energy
        E -= 0.25 * Tensor<0>(ovov.contract(t, Eigen::array<Eigen::IndexPair<int>, 4>{a1, a2, a3, a4}))(0);

        // print the MP2 timing
        std::printf("MP2 CORRELATION ENERGY TIMING: %s\n", Timer::Format(Timer::Elapsed(tp)).c_str());
    }

    // add the MP3 correlation energy
    if (tp = Timer::Now(); program.get<int>("-o") > 2) {
        // define the contraction axes
        Eigen::IndexPair<int> a11(1, 0), a12(3, 2), a21(0, 0), a22(1, 2), a23(2, 1), a24(3, 3);
        Eigen::IndexPair<int> b11(0, 1), b12(2, 3), b21(0, 1), b22(1, 3), b23(2, 0), b24(3, 2);
        Eigen::IndexPair<int> c11(0, 1), c12(3, 2), c21(0, 1), c22(1, 2), c23(2, 0), c24(3, 3);

        // add the MP3 correlation energy
        E += 0.125 * Tensor<0>(t.contract(vvvv, Eigen::array<Eigen::IndexPair<int>, 2>{a11, a12}).contract(t, Eigen::array<Eigen::IndexPair<int>, 4>{a21, a22, a23, a24}))(0);
        E += 0.125 * Tensor<0>(t.contract(oooo, Eigen::array<Eigen::IndexPair<int>, 2>{b11, b12}).contract(t, Eigen::array<Eigen::IndexPair<int>, 4>{b21, b22, b23, b24}))(0);
        E -=         Tensor<0>(t.contract(oovv, Eigen::array<Eigen::IndexPair<int>, 2>{c11, c12}).contract(t, Eigen::array<Eigen::IndexPair<int>, 4>{c21, c22, c23, c24}))(0);

        // print the MP2 timing
        std::printf("MP3 CORRELATION ENERGY TIMING: %s\n", Timer::Format(Timer::Elapsed(tp)).c_str());
    }

    // print the final energy
    std::printf("\nFINAL SINGLE POINT ENERGY: %.14f\n\nTOTAL TIME: %s\n", E + system.nuclearRepulsion(), Timer::Format(Timer::Elapsed(start)).c_str());
}
