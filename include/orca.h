#pragma once

#include "method.h"

class Orca : public Method<Orca> {
    friend class Method<Orca>;
public:
    struct Options {Options(){};
        // dynamics structure
        struct Dynamics {
            struct Berendsen {bool enabled=false; double tau = 1, temp=298.15, timeout=10;} berendsen={};
            int iters=100; double step=1.0;
        } dynamics={};

        // variables
        std::filesystem::path interface; std::string method;
    };
public:
    // constructors and destructors
    Orca(const Options& opt = {}) : opt(opt) {}

    // overriden virtual derivatives
    Result gradient(const System& system, const Integrals& ints = {}, Result res = {}, bool print = true) const override;
    Result run(const System& system, const Integrals& ints = {}, Result res = {}, bool print = true) const override;

private:
    Options opt;
};
