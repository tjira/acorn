#include "input.h"

nlohmann::json default_input = R"(
{
    "system" : {
        "path" : "molecule.xyz",
        "basis" : "sto-3g",
        "charge" : 0,
        "multiplicity" : 1
    },
    "wavefunction" : {
        "guess" : ["-(x-1)^2"],
        "dimension" : 1,
        "grid_limits" : [-8, 8],
        "grid_points" : 8,
        "mass" : 1.0,
        "momentum" : 0.0
    },
    "integral" : {
        "precision" : 1e-12,
        "data_export" : {
            "hamiltonian" : false,
            "coulomb" : false,
            "overlap" : false
        }
    },
    "hartree_fock" : {
        "generalized" : false,
        "diis_size" : 5,
        "max_iter" : 100,
        "threshold" : 1e-8,
        "moller_plesset" : {
            "order" : 2
        },
        "configuration_interaction" : {
            "excitation" : []
        },
        "coupled_cluster" : {
            "excitation" : [1, 2],
            "max_iter" : 100,
            "threshold" : 1e-8,
            "perturbation" : [3]
        }
    },
    "quantum_dynamics" : {
        "potential" : [["0.5*x^2"]],
        "imaginary" : 0,
        "real" : 0,
        "iterations" : 1000,
        "time_step" : 1.0,
        "data_export" : {
            "diabatic_wavefunction" : false,
            "adiabatic_wavefunction" : false,
            "diabatic_density" : false,
            "adiabatic_density" : false,
            "diabatic_population" : false,
            "adiabatic_population" : false,
            "diabatic_potential" : false,
            "adiabatic_potential" : false,
            "energy" : false,
            "position" : false,
            "momentum" : false,
            "acf" : false
        }
    },
    "classical_dynamics" : {
        "potential" : [["0.5*x^2"]],
        "iterations" : 1000,
        "seed" : 1,
        "adiabatic" : false,
        "trajectories" : 100,
        "time_step" : 1.0,
        "log_interval" : 1,
        "data_export" : {
            "diabatic_population" : false,
            "adiabatic_population" : false,
            "energy" : false,
            "position" : false,
            "momentum" : false
        }
    }
}
)"_json;

int nthread = 1;
