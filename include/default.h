#pragma once

#include <nlohmann/json.hpp>

inline nlohmann::json rhfopt = {
    {"maxiter", 100}, {"thresh", 1e-8}
};
