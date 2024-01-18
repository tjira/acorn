#pragma once

#include "integral.h"
#include "timer.h"

struct Result {
    struct RestrictedHartreeFock {
        double E; Vector<> eps; Matrix<> C, D, G, H;
    } rhf;
    struct RestrictedMollerPlesset {
        double Ecorr; Matrix<> G, H;
    } rmp;
    double Etot; Matrix<> G, H;
};

template <class M>
class Method {
public:
    // general methods
    void dynamics(System system, const Integrals& ints, Result res, bool print = true) const;

    // static methods
    static Vector<> frequency(const System& system, const Matrix<>& H);

    // virtual derivatives
    virtual Result gradient(const System& system, const Integrals& ints, Result res, bool print = true) const;
    virtual Result hessian(const System& system, const Integrals& ints, Result res, bool print = true) const;

    // virtual method functions
    virtual Result run(const System&, const Integrals&, Result res, bool) const {return res;};

    // child type converter
    const M* get() const {return static_cast<const M*>(this);}

private:
    virtual Result run(const System&, Result res, bool) const {return res;};
};
