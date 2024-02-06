#pragma once

#include "method.h"

class Orca : public Method<Orca> {
    friend class Method<Orca>;
public:
    struct Options {Options(){};
        // structures
        struct Dynamics {int iters=100; double step=1.0; std::string folder;} dynamics={};

        // variables
        std::string method; std::filesystem::path interface, folder;
    };
public:
    // constructors and destructors
    Orca(const Options& opt = {}) : opt(opt) {}

    // overriden virtual derivatives
    Result gradient(const System& system, const Integrals& ints, Result res, bool print = true) const override;

private:
    Options opt;
};
