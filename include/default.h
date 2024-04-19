#pragma once

#include <string>

inline std::string rcioptstr = R"({
    "excitations" : [], "nstate" : 3, "state" : 1,
    "dynamics" : {"iters" : 100, "step" : 20, "berendsen" : {"tau" : 15, "temp" : 0, "timeout" : 20}}, "gradient" : {"step" : 1e-5}, "hessian" : {"step" : 1e-5},
    "export" : {"coulombms" : false, "kineticms" : false, "nuclearms" : false, "hcorems" : false, "hamiltonian" : false, "energies" : false},
    "print" : {"coulombms" : false, "kineticms" : false, "nuclearms" : false, "hcorems" : false, "hamiltonian" : false, "energies" : true}
})";

inline std::string intoptstr = R"({
    "export" : {"kinetic" : false, "nuclear" : false, "overlap" : false, "coulomb" : false, "dkinetic" : false, "dnuclear" : false, "doverlap" : false, "dcoulomb" : false},
    "print" : {"kinetic" : false, "nuclear" : false, "overlap" : false, "coulomb" : false, "dkinetic" : false, "dnuclear" : false, "doverlap" : false, "dcoulomb" : false}
})";

inline std::string mdloptstr = R"({
    "mass" : 1, "ngrid" : 512, "limits" : [-16, 16], "variables" : ["x"]
})";

inline std::string moloptstr = R"({
    "charge" : 0, "multiplicity" : 1
})";

inline std::string msaoptstr = R"({
    "real" : false, "step" : 0.1, "guess" : "-x^2", "nstate" : 1, "iters" : 1000, "savewfn" : false,
    "optimize" : {"step" : 1, "iters" : 1000},
    "spectrum" : {"normalize" : false, "zpesub" : false, "potential" : "", "window" : "1", "zeropad" : 0},
    "dynamics" : {"iters" : 100, "step" : 20, "state" : 1, "position" : [0], "velocity" : [0], "gradient" : []}
})";

inline std::string msnoptstr = R"({
    "step" : 0.1, "guess" : ["-x^2"], "iters" : 1000, "savewfn" : false, "cap" : "0", "momentum" : 0, "adiabatic" : true,
    "dynamics" : {"iters" : 100, "step" : 20, "state" : 1, "position" : [0], "velocity" : [0], "gradient" : []}
})";

inline std::string rmpoptstr = R"({
    "order" : 2,
    "dynamics" : {"iters" : 100, "step" : 20, "berendsen" : {"tau" : 15, "temp" : 0, "timeout" : 20}}, "gradient" : {"step" : 1e-5}, "hessian" : {"step" : 1e-5},
    "export" : {"coulombmo" : false},
    "print" : {"coulombmo" : false}
})";

inline std::string rhfoptstr = R"({
    "maxiter" : 1000, "thresh" : 1e-12, "mulliken" : false,
    "dynamics" : {"iters" : 100, "step" : 20, "berendsen" : {"tau" : 15, "temp" : 0, "timeout" : 20}}, "gradient" : {"step" : 1e-5}, "hessian" : {"step" : 1e-5},
    "export" : {"coef" : false, "density" : false, "orben" : false, "hcore" : false},
    "print" : {"coef" : false, "density" : false, "orben" : false, "hcore" : false}
})";

inline std::string uhfoptstr = R"({
    "maxiter" : 1000, "thresh" : 1e-12,
    "dynamics" : {"iters" : 100, "step" : 20, "berendsen" : {"tau" : 15, "temp" : 0, "timeout" : 20}}, "gradient" : {"step" : 1e-5}, "hessian" : {"step" : 1e-5},
    "export" : {"coefa" : false, "coefb" : false, "densitya" : false, "densityb" : false, "orbena" : false, "orbenb" : false, "hcore" : false},
    "print" : {"coefa" : false, "coefb" : false, "densitya" : false, "densityb" : false, "orbena" : false, "orbenb" : false, "hcore" : false}
})";

inline std::string umpoptstr = R"({
    "order" : 2,
    "dynamics" : {"iters" : 100, "step" : 20, "berendsen" : {"tau" : 15, "temp" : 0, "timeout" : 20}}, "gradient" : {"step" : 1e-5}, "hessian" : {"step" : 1e-5},
    "export" : {"coulombmo" : false},
    "print" : {"coulombmo" : false}
})";

inline std::string bgloptstr = R"({
    "interface" : "bagel.sh", "nstate" : 3, "state" : 1,
    "dynamics" : {"iters" : 100, "step" : 20, "berendsen" : {"tau" : 15, "temp" : 0, "timeout" : 20}}
})";

inline std::string orcoptstr = R"({
    "interface" : "orca.sh",
    "dynamics" : {"iters" : 100, "step" : 20, "berendsen" : {"tau" : 15, "temp" : 0, "timeout" : 20}}
})";
