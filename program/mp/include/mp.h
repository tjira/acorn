#pragma once

#define OCC "abcdefghijklmnopqrstuvwxyz"
#define VRT "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

#include "tensor.h"

namespace Acorn {
    namespace MBPT {
        double evaluate(const std::string& contrstr, const torch::Tensor& Jmsa, const torch::Tensor& Ems, int nos, int o);
        struct Options {
            int order;
        };
        inline int nthread = 1;
    }
}

#define SPLIT(S, D) [](const std::string& s, char d) {std::vector<std::string> r; std::stringstream ss(s); std::string l; while (getline(ss, l, d)) {r.push_back(l);} return r;}(S, D)
