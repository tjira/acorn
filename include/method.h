#pragma once

#include "integral.h"
#include "result.h"
#include "timer.h"

template <class M>
class Method {
public:
    // general methods
    void dynamics(System system, Integrals ints, Result res, bool print = true) const;

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
    virtual std::tuple<Result, Integrals> run(const System&, Result res, bool) const {return {res, {}};};
};
