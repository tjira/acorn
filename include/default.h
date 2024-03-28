#pragma once

#include <string>

inline std::string rcioptstr = R"({
    "excitations" : [],
    "dynamics" : {"iters" : 100, "step" : 20, "berendsen" : {"tau" : 15, "temp" : 0, "timeout" : 20}}, "gradient" : {"step" : 1e-5}, "hessian" : {"step" : 1e-5},
    "export" : {"coulombms" : false, "kineticms" : false, "nuclearms" : false, "hcorems" : false, "hamiltonian" : false, "energies" : false},
    "print" : {"coulombms" : false, "kineticms" : false, "nuclearms" : false, "hcorems" : false, "hamiltonian" : false, "energies" : false}
})";

inline std::string intoptstr = R"({
    "export" : {"kinetic" : false, "nuclear" : false, "overlap" : false, "coulomb" : false, "dkinetic" : false, "dnuclear" : false, "doverlap" : false, "dcoulomb" : false},
    "print" : {"kinetic" : false, "nuclear" : false, "overlap" : false, "coulomb" : false, "dkinetic" : false, "dnuclear" : false, "doverlap" : false, "dcoulomb" : false}
})";

inline std::string mdloptstr = R"({
    "mass" : 1, "ngrid" : 512, "limits" : [[-16, 16]]
})";

inline std::string moloptstr = R"({
    "charge" : 0, "multiplicity" : 1
})";

inline std::string msaoptstr = R"({
    "real" : false, "step" : 0.1, "optimize" : false, "guess" : "-x^2", "nstate" : 1, "iters" : 1000, "savewfn" : false,
    "spectrum" : {"normalize" : false, "zpesub" : false, "potential" : "", "window" : ""},
    "dynamics" : {"iters" : 100, "step" : 20, "berendsen" : {"tau" : 15, "temp" : 0, "timeout" : 20}}
})";

inline std::string msnoptstr = R"({
    "step" : 0.1, "guess" : ["-x^2"], "iters" : 1000, "savewfn" : false,
    "dynamics" : {"iters" : 100, "step" : 20, "berendsen" : {"tau" : 15, "temp" : 0, "timeout" : 20}}
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
    "interface" : "orca.sh",
    "dynamics" : {"iters" : 100, "step" : 20, "berendsen" : {"tau" : 15, "temp" : 0, "timeout" : 20}}
})";

inline std::string orcoptstr = R"({
    "interface" : "orca.sh",
    "dynamics" : {"iters" : 100, "step" : 20, "berendsen" : {"tau" : 15, "temp" : 0, "timeout" : 20}}
})";
