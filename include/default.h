#pragma once

#include <nlohmann/json.hpp>

#define CFREQ 5140.486777894163
#define BOHR2A 0.529177210903
#define A2BOHR 1 / BOHR2A
#define BOLTZMANN 3.166811429e-6
#define AU2FS 0.02418884254

inline int nthread;

inline nlohmann::json mpopt = {
    {"order", 2}, {"gradient", false}
};

inline nlohmann::json rhfopt = {
    {"maxiter", 100}, {"thresh", 1e-8}, {"gradient", false}
};
