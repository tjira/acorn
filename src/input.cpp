#include "input.h"

nlohmann::json default_input = R"(
{
    "system" : {
        "path" : "molecule.xyz",
        "basis" : "sto-3g"
    },
    "wavefunction" : {
        "guess" : ["-(x-1)^2"],
        "dimension" : 1,
        "grid_limits" : [-8, 8],
        "grid_points" : 8,
        "mass" : 1.0,
        "momentum" : 0.0
    },
    "model" : {
        "potential" : [["0.5*x^2"]],
        "dimension" : 1,
        "mass" : 1.0,
        "momentum" : 0.0
    },
    "integral" : {
        "precision" : 1e-12
    },
    "hartree_fock" : {
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
        "align_wavefunction" : true,
        "data_export" : {
            "diabatic_wavefunction" : false,
            "adiabatic_wavefunction" : false,
            "diabatic_density" : false,
            "adiabatic_density" : false,
            "diabatic_potential" : false,
            "adiabatic_potential" : false,
            "energy" : false,
            "position" : false,
            "momentum" : false,
            "acf" : false
        }
    }
}
)"_json;
