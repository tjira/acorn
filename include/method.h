#pragma once

#include "integral.h"
#include "timer.h"

struct Result {
    struct RestrictedHartreeFock {
        double E; Vector<> eps; Matrix<> C, D;
    } rhf;
    struct RestrictedMollerPlesset {
        double Ecorr;
    } rmp;
    double Etot; Matrix<> G, H;
};

template <class M>
class Method {
public:
    // general methods
    Result gradient(const System& system, const Integrals& ints, Result res, bool print = true) const;
    Result hessian(const System& system, const Integrals& ints, Result res, bool print = true) const;
    void dynamics(System system, const Integrals& ints, Result res, bool print = true) const;

    // static methods
    static Vector<> frequency(const System& system, const Matrix<>& H);

    // virtual method functions
    virtual Result run(const System&, const Integrals&, Result, bool) const = 0;
    virtual Result run(const System&, Result, bool) const = 0;

    // child type converter
    const M* get() const {return static_cast<const M*>(this);}
};

#include <nlohmann/json.hpp>
