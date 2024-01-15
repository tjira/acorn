#pragma once

#include <nlohmann/json.hpp>

#define CFREQ 5140.486777894163
#define BOHR2A 0.529177210903
#define A2BOHR 1 / BOHR2A
#define BOLTZMANN 3.166811429e-6
#define AU2FS 0.02418884254

inline int nthread;

inline std::unordered_map<int, double> masses = {
    {1, 01.007840},
    {6, 12.011000},
    {7, 14.006700},
    {8, 15.999000},
    {9, 18.998403},
};

inline nlohmann::json rmpopt = {
    {"order", 2}, {"dynamics", {{"iters", 100}, {"output", "trajectory.xyz"}, {"step", 1}}}, {"gradient", {{"step", 1e-5}}}, {"hessian", {{"step", 1e-5}}}
};

inline nlohmann::json rhfopt = {
    {"maxiter", 1000}, {"thresh", 1e-12}, {"dynamics", {{"iters", 100}, {"output", "trajectory.xyz"}, {"step", 1}}}, {"gradient", {{"step", 1e-5}}}, {"hessian", {{"step", 1e-5}}}
};
