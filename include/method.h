#pragma once

#include "integral.h"

class Method {
public:
    struct Options{
        struct Gradient {double step;} gradient; struct Hessian {double step;} hessian;
        struct Dynamics {int iters; double step; std::string output;} dynamics;
    } opt;
    struct Result {
        struct RestrictedHartreeFock {
            double E; Vector<> eps; Matrix<> C, D;
        } rhf;
        struct RestrictedMollerPlesset {
            double Ecorr;
        } rmp;
        double Etot; Matrix<> G, H;
    };
public:
    // constructors and destructors
    Method(const Method::Options& opt) : opt(opt) {}

    // general methods
    Result gradient(const System& system, const Integrals& ints, Result res, bool print = true) const;
    Result hessian(const System& system, const Integrals& ints, Result res, bool print = true) const;
    void dynamics(System system, Integrals ints, Result res, bool print = true) const;

    // static methods
    static Vector<> frequency(const System& system, const Matrix<>& H);

    // virtual functions
    virtual double energy(const System&, Result) const = 0;
    virtual Result run(const System&, const Integrals&, Result, bool) const = 0;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Method::Options::Dynamics, iters, step, output);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Method::Options::Gradient, step);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Method::Options::Hessian, step);
