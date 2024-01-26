#pragma once

#include <string>

inline std::string intoptstr = R"({
    "export" : {"kinetic" : false, "nuclear" : false, "overlap" : false, "coulomb" : false},
    "print" : {"kinetic" : false, "nuclear" : false, "overlap" : false, "coulomb" : false}
})";

inline std::string mdloptstr = R"({
    "mass" : 1, "ngrid" : 512, "limits" : [[-16, 16]]
})";

inline std::string msaoptstr = R"({
    "real" : false, "step" : 0.1, "optimize" : true, "guess" : "-x^2", "nstate" : 3, "iters" : 1000, "thresh" : 1e-12, "dynamics" : {"iters" : 100, "output" : "trajectory.xyz", "step" : 1}
})";

inline std::string msnoptstr = R"({
    "step" : 0.1, "guess" : ["-x^2"], "iters" : 1000, "dynamics" : {"iters" : 100, "output" : "trajectory.xyz", "step" : 1}
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
