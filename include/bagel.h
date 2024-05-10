#pragma once

#include "method.h"

class Bagel : public Method<Bagel> {
    friend class Method<Bagel>;
public:
    struct Options {Options(){};
        // dynamics structure
        struct Dynamics {
            int iters=100; double step=1.0;
        } dynamics={};

        // variables
        std::filesystem::path interface; std::string method; int nstate=2, state=1;
    };
public:
    // constructors and destructors
    Bagel(const Options& opt = {}) : opt(opt) {}

    // overriden virtual derivatives
    Result gradient(const System& system, const Integrals& ints = {}, Result res = {}, bool print = true) const override;
    Result run(const System& system, const Integrals& ints = {}, Result res = {}, bool print = true) const override;

private:
    Options opt;
};
