#pragma once

#include <string>

inline std::string intoptstr = R"({
    "export" : {"kinetic" : false, "nuclear" : false, "overlap" : false, "coulomb" : false},
    "print" : {"kinetic" : false, "nuclear" : false, "overlap" : false, "coulomb" : false}
})";

inline std::string msvoptstr = R"({
    "real" : false, "step" : 0.1, "optimize" : true, "nstate" : 3, "iters" : 1000, "thresh" : 1e-12, "dynamics" : {"iters" : 100, "output" : "trajectory.xyz", "step" : 1}
})";

inline std::string rmpoptstr = R"({
    "order" : 2, "dynamics" : {"iters" : 100, "output" : "trajectory.xyz", "step" : 1}, "gradient" : {"step" : 1e-5}, "hessian" : {"step" : 1e-5}
})";

inline std::string rhfoptstr = R"({
    "maxiter" : 1000, "thresh" : 1e-12, "dynamics" : {"iters" : 100, "output" : "trajectory.xyz", "step" : 1}, "gradient" : {"step" : 1e-5}, "hessian" : {"step" : 1e-5}
})";

inline std::string orcoptstr = R"({
    "interface" : "orca.sh", "method" : "hf", "dynamics" : {"iters" : 100, "output" : "trajectory.xyz", "step" : 1}
})";
