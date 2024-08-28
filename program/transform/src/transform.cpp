#include "transform.h"

Tensor<4> Acorn::Transform::CoulombSpatial(const Tensor<4>& Jao, Matrix& Cmo) {
    // perform the transform
    Tensor<4> J01 = TENSORMAP(Cmo).contract(Jao, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 0)});
    Tensor<4> J02 = TENSORMAP(Cmo).contract(J01, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 1)});
    Tensor<4> J03 = TENSORMAP(Cmo).contract(J02, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 2)});
    Tensor<4> Jmo = TENSORMAP(Cmo).contract(J03, Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(0, 3)});

    // return the Coulomb integrals in molecular orbital basis
    return Jmo;
}

Tensor<4> Acorn::Transform::CoulombSpin(const Tensor<4>& Jao, const Matrix& Cmo) {
    // create the tiling matrix P that repeats the MO columns 2 times and define the spin mask
    Matrix P = Matrix::NullaryExpr(Cmo.rows(), 2 * Cmo.rows(), std::function<double(int, int)>([](int i, int j) {return j == 2 * i || j == 2 * i + 1;})), MN(2 * Cmo.cols(), 2 * Cmo.cols());

    // initialize the spin mask
    MN << Matrix::NullaryExpr(Cmo.rows(), 2 * Cmo.cols(), std::function<double(int, int)>([](int, int j) {return 1 - j % 2;})),
          Matrix::NullaryExpr(Cmo.rows(), 2 * Cmo.cols(), std::function<double(int, int)>([](int, int j) {return 0 + j % 2;}));

    // transform the wfn coefficients to the spin basis
    Matrix Cms = (Cmo * P).replicate<2, 1>().cwiseProduct(MN);

    // return the transformed matrix
    return CoulombSpatial(Eigen::Kron(Matrix::Identity(2, 2), Eigen::Kron(Matrix::Identity(2, 2), Jao).shuffle(Eigen::array<int, 4>{3, 2, 1, 0})), Cms);
}

Matrix Acorn::Transform::SingleSpatial(const Matrix& Aao, Matrix& Cmo) {return Cmo.transpose() * Aao * Cmo;}

Matrix Acorn::Transform::SingleSpin(const Matrix& Aao, const Matrix& Cmo) {
    // create the tiling matrix P that repeats the MO columns 2 times and define the spin mask
    Matrix P = Matrix::NullaryExpr(Cmo.rows(), 2 * Cmo.rows(), std::function<double(int, int)>([](int i, int j) {return j == 2 * i || j == 2 * i + 1;})), MN(2 * Cmo.cols(), 2 * Cmo.cols());

    // initialize the spin mask
    MN << Matrix::NullaryExpr(Cmo.rows(), 2 * Cmo.cols(), std::function<double(int, int)>([](int, int j) {return 1 - j % 2;})),
          Matrix::NullaryExpr(Cmo.rows(), 2 * Cmo.cols(), std::function<double(int, int)>([](int, int j) {return 0 + j % 2;}));

    // transform the wfn coefficients to the spin basis
    Matrix Cms = (Cmo * P).replicate<2, 1>().cwiseProduct(MN);

    // return the transformed matrix
    return SingleSpatial(Eigen::kroneckerProduct(Matrix::Identity(2, 2), Aao), Cms);
}

void Acorn::Transform::run(const Options& opt, std::vector<timepoint>& timers) {
    // define all the matrices and tensors used throughout the program
    Matrix C, E, V, T, S, F, Vmo, Tmo, Smo, Fmo, Ems, Vms, Tms, Sms, Fms; Tensor<4> J, Jmo, Jms;

    // load the coefficient matrix and integrals in AO basis from disk
    if (opt.orbital || opt.spinorbital) {

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
        Ems = Vector::NullaryExpr(2 * E.rows(), [&](int i){return E(i / 2);});

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
        Eigen::Write("V_MO.mat", Vmo);
        Eigen::Write("T_MO.mat", Tmo);
        Eigen::Write("S_MO.mat", Smo);
        Eigen::Write("J_MO.mat", Jmo);
        Eigen::Write("F_MO.mat", Fmo);

        // print the time for integral writing
    }

    // integrals in molecular spatial orbital basis writing
    if (opt.spinorbital) {

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
        std::cout << eltime(timers.at(1)) << std::endl;
    }

    // print the total time
    std::cout << (opt.orbital || opt.spinorbital ? "\n" : "") << "TOTAL TIME: " << eltime(timers.at(0)) << std::endl;
}
