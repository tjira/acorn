#include "transform.h"

torch::Tensor Acorn::Transform::CoulombSpatial(const torch::Tensor& Jao, torch::Tensor& Cmo) {return torch::einsum("ip,jq,ijkl,kr,ls->pqrs", {Cmo, Cmo, Jao, Cmo, Cmo});}

torch::Tensor Acorn::Transform::CoulombSpin(const torch::Tensor& Jao, const torch::Tensor& Cmo) {
    // create the tile matrix P
    torch::Tensor P = torch::eye(Cmo.sizes().at(0), Cmo.sizes().at(1), torch::kDouble).repeat_interleave(2, 1);

    // create the spin mask matrices
    torch::Tensor N = torch::arange(2, torch::kDouble).repeat({Cmo.sizes().at(0), Cmo.sizes().at(1)}); torch::Tensor M = (N + 1) % 2;

    // transform the wfn coefficients to the spin basis
    torch::Tensor Cms = torch::cat({Cmo.mm(P), Cmo.mm(P)}) * torch::cat({M, N});

    // return the transformed matrix
    return CoulombSpatial(torch::kron(torch::eye(2, 2, torch::kDouble), torch::kron(torch::eye(2, 2, torch::kDouble), Jao).swapaxes(0, 3).swapaxes(1, 2).contiguous()), Cms);
}

torch::Tensor Acorn::Transform::SingleSpatial(const torch::Tensor& Aao, torch::Tensor& Cmo) {return Cmo.t().mm(Aao).mm(Cmo);}

torch::Tensor Acorn::Transform::SingleSpin(const torch::Tensor& Aao, const torch::Tensor& Cmo) {
    // create the tile matrix P
    torch::Tensor P = torch::eye(Cmo.sizes().at(0), Cmo.sizes().at(1), torch::kDouble).repeat_interleave(2, 1);

    // create the spin mask matrices
    torch::Tensor N = torch::arange(2, torch::kDouble).repeat({Cmo.sizes().at(0), Cmo.sizes().at(1)}); torch::Tensor M = (N + 1) % 2;

    // transform the wfn coefficients to the spin basis
    torch::Tensor Cms = torch::cat({Cmo.mm(P), Cmo.mm(P)}) * torch::cat({M, N});

    // return the transformed matrix
    return SingleSpatial(torch::kron(torch::eye(2, 2, torch::kDouble), Aao), Cms);
}

void Acorn::Transform::run(const Options& opt, std::vector<timepoint>& timers) {
    // define all the matrices and tensors used throughout the program
    torch::Tensor C, E, F, J, T, V, S, Fmo, Jmo, Tmo, Vmo, Smo, Ems, Fms, Jms, Tms, Vms, Sms;

    // load the coefficient matrix and integrals in AO basis from disk
    if (opt.orbital || opt.spinorbital) {

        // start the timer for integral loading
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // print the header of integral loading
        std::cout << "SYSTEM AND INTEGRALS IN AO BASIS READING: " << std::flush;

        // load the integrals in AO basis
        E = torch::ReadTensor("E_MO.mat");
        C = torch::ReadTensor("C_MO.mat");
        V = torch::ReadTensor("V_AO.mat");
        T = torch::ReadTensor("T_AO.mat");
        S = torch::ReadTensor("S_AO.mat");
        J = torch::ReadTensor("J_AO.mat");
        F = torch::ReadTensor("F_AO.mat");

        // print the time for integral loading
        std::cout << eltime(timers.at(1)) << std::endl;
    }

    // print a new line
    if (opt.orbital || opt.spinorbital) std::cout << std::endl;

    // integrals in MO basis calculation
    if (opt.orbital) {
        
        // start the timer for integral transformation
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // print the header of integral transformation
        std::cout << "INTEGRAL TRANSFORMS TO MO BASIS: " << std::flush;

        // calculate the integrals in MO basis
        Jmo = Acorn::Transform::CoulombSpatial(J, C);
        Vmo = Acorn::Transform::SingleSpatial (V, C);
        Tmo = Acorn::Transform::SingleSpatial (T, C);
        Smo = Acorn::Transform::SingleSpatial (S, C);
        Fmo = Acorn::Transform::SingleSpatial (F, C);

        // print the time for integral transformation
        std::cout << eltime(timers.at(1)) << std::endl;
    }

    // integrals in MO basis calculation
    if (opt.spinorbital) {
    
        // start the timer for integral transformation
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // print the header of integral transformation
        std::cout << "INTEGRAL TRANSFORMS TO MS BASIS: " << std::flush;

        // calculate the integrals in MS basis
        Jms = Acorn::Transform::CoulombSpin(J, C);
        Vms = Acorn::Transform::SingleSpin (V, C);
        Tms = Acorn::Transform::SingleSpin (T, C);
        Sms = Acorn::Transform::SingleSpin (S, C);
        Fms = Acorn::Transform::SingleSpin (F, C);
        Ems = E.repeat_interleave(2);

        // print the time for integral transformation
        std::cout << eltime(timers.at(1)) << std::endl;
    }

    // integrals in molecular spatial orbital basis writing
    if (opt.orbital) {

        // start the timer for integral writing
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // print the header of integral writing
        std::cout << "INTEGRALS IN MO BASIS WRITING:   " << std::flush;

        // save the integrals to disk
        torch::WriteTensor("V_MO.mat", Vmo);
        torch::WriteTensor("T_MO.mat", Tmo);
        torch::WriteTensor("S_MO.mat", Smo);
        torch::WriteTensor("J_MO.mat", Jmo);
        torch::WriteTensor("F_MO.mat", Fmo);

        // print the time for integral writing
        std::cout << eltime(timers.at(1)) << std::endl;
    }

    // integrals in molecular spatial orbital basis writing
    if (opt.spinorbital) {

        // start the timer for integral writing
        timers.at(1) = std::chrono::high_resolution_clock().now();

        // print the header of integral writing
        std::cout << "INTEGRALS IN MS BASIS WRITING:   " << std::flush;

        // save the integrals to disk
        torch::WriteTensor("V_MS.mat", Vms);
        torch::WriteTensor("T_MS.mat", Tms);
        torch::WriteTensor("S_MS.mat", Sms);
        torch::WriteTensor("J_MS.mat", Jms);
        torch::WriteTensor("E_MS.mat", Ems);
        torch::WriteTensor("F_MS.mat", Fms);

        // print the time for integral writing
        std::cout << eltime(timers.at(1)) << std::endl;
    }

    // print the total time
    std::cout << (opt.orbital || opt.spinorbital ? "\n" : "") << "TOTAL TIME: " << eltime(timers.at(0)) << std::endl;
}
