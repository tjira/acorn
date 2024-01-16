#pragma once

#include <nlohmann/json.hpp>

inline nlohmann::json intopt = {
    {"export", {{"kinetic", false}, {"nuclear", false}, {"overlap", false}, {"coulomb", false}}},
    {"print", {{"kinetic", false}, {"nuclear", false}, {"overlap", false}, {"coulomb", false}}}
};

inline nlohmann::json rmpopt = {
    {"order", 2}, {"dynamics", {{"iters", 100}, {"output", "trajectory.xyz"}, {"step", 1}}}, {"gradient", {{"step", 1e-5}}}, {"hessian", {{"step", 1e-5}}}
};

inline nlohmann::json rhfopt = {
    {"maxiter", 1000}, {"thresh", 1e-12}, {"dynamics", {{"iters", 100}, {"output", "trajectory.xyz"}, {"step", 1}}}, {"gradient", {{"step", 1e-5}}}, {"hessian", {{"step", 1e-5}}}
};
