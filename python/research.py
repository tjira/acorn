#!/usr/bin/env python

import copy as cp, glob as gl, shutil as sh, json as js, numpy as np, os

dollar = "$"  if os.name == "nt" else  "\\$"
esc    = "\\" if os.name == "nt" else "\\\\"

KTSH_COUPLING = 1; KTSH_HO = 1; SH_COMPARE = 1; LZ_COMPARE = 1; URACIL_LVC = 1

TRAJS, NGRID = 10000, 8192

# KTSH TDC TESTING ON MODEL POTENTIALS ======================================================================================================================================================

if KTSH_COUPLING:

    HAMILTONIAN = {
        "P1" : {
            "dims" : 1, "states" : 2,
            "matrix" : [
                ["sgn(r1) * 0.01 * (1 - exp(-(1.6 * sgn(r1) * r1)))", "0.005 * exp(-(r1^2))"                              ],
                ["0.005 * exp(-(r1^2))",                              "-sgn(r1) * 0.01 * (1 - exp(-(1.6 * sgn(r1) * r1)))"],
            ]
        },
        "P2" : {
            "dims" : 1, "states" : 2,
            "matrix" : [
                ["sgn(r1) * 0.01 * (1 - exp(-(1.6 * sgn(r1) * r1)))", "0.0005 / (200 * r1^2 + 1)"                         ],
                ["0.0005 / (200 * r1^2 + 1)",                         "-sgn(r1) * 0.01 * (1 - exp(-(1.6 * sgn(r1) * r1)))"],
            ]
        },
        "P3" : {
            "dims" : 1, "states" : 3,
            "matrix" : [
                ["sgn(r1) * 0.01 * (1 - exp(-(1.6 * sgn(r1) * r1)))", "0.003 * exp(-(r1^2))", "0"                                                 ],
                ["0.003 * exp(-(r1^2))",                              "0",                    "0.003 * exp(-(r1^2))"                              ],
                ["0",                                                 "0.003 * exp(-(r1^2))", "-sgn(r1) * 0.01 * (1 - exp(-(1.6 * sgn(r1) * r1)))"],
            ]
        },
        "P4" : {
            "dims" : 1, "states" : 3,
            "matrix" : [
                ["sgn(r1) * 0.01 * (1 - exp(-(1.6 * sgn(r1) * r1)))", "0.0005 / (100 * r1^2 + 1)", "0"                                                 ],
                ["0.0005 / (100 * r1^2 + 1)",                         "0",                         "0.0005 / (100 * r1^2 + 1)"                         ],
                ["0",                                                 "0.0005 / (100 * r1^2 + 1)", "-sgn(r1) * 0.01 * (1 - exp(-(1.6 * sgn(r1) * r1)))"],
            ]
        }
    }

    IS = {
        "P1" : [1, 1],
        "P2" : [1, 1],
        "P3" : [2, 2],
        "P4" : [2, 2],
    }

    NS = {
        "P1" : 2,
        "P2" : 2,
        "P3" : 3,
        "P4" : 3,
    }

    TDCLIM = {
        "P1" : [920, 1520],
        "P2" : [920, 1520],
        "P3" : [920, 1520],
        "P4" : [920, 1520],
    }

    DS, TS = ["P1", "P2"], ["P3", "P4"]

    for potential in DS + TS:

        dp_input = {"model_potential" : {"adiabatic" : False, "limits" : [-16, 16], "output" : f"KTSH_COUPLING_POTENTIAL_DIABATIC_{potential}.mat",  "points" : 8192, "hamiltonian" : HAMILTONIAN[potential]}}
        ap_input = {"model_potential" : {"adiabatic" : True,  "limits" : [-16, 16], "output" : f"KTSH_COUPLING_POTENTIAL_ADIABATIC_{potential}.mat", "points" : 8192, "hamiltonian" : HAMILTONIAN[potential]}}

        qd_input = {
            "quantum_dynamics" : {
                "log_intervals" : {
                    "iteration" : 50
                },

                "grid" : {
                    "limits" : [-24, 48], "points" : NGRID
                },

                "hamiltonian" : HAMILTONIAN[potential],

                "initial_conditions" : {
                    "position" : [-10], "momentum" : [15], "gamma" : 2, "state" : IS[potential][0], "mass" : 2000
                },

                "write" : {
                    "population" : f"KTSH_COUPLING_POPULATION_{potential}_EXACT.mat"
                },

                "adiabatic" : True,
                "iterations" : 300,
                "mode" : [0, 1],
                "time_step" : 10
            }
        }

        cd_input = {
            "classical_dynamics" : {
                "log_intervals" : {
                    "iteration" : 1000, "trajectory" : 1000
                },

                "hamiltonian" : HAMILTONIAN[potential],

                "initial_conditions" : {
                    "position_mean" : qd_input["quantum_dynamics"]["initial_conditions"]["position"], "position_std"  : [0.5],
                    "momentum_mean" : qd_input["quantum_dynamics"]["initial_conditions"]["momentum"], "momentum_std"  : [1.0],
                    "state" : [int(i == IS[potential][1]) for i in range(NS[potential])],
                    "mass" : [qd_input["quantum_dynamics"]["initial_conditions"]["mass"]]
                },

                "write" : {},

                "adiabatic" : qd_input["quantum_dynamics"]["adiabatic"],
                "iterations" : qd_input["quantum_dynamics"]["iterations"] * 10,
                "time_step" : qd_input["quantum_dynamics"]["time_step"] / 10,
                "trajectories" : TRAJS
            }
        }

        fs_input = cp.deepcopy(cd_input); fs_input["classical_dynamics"]["fewest_switches"] = {}
        kt_input = cp.deepcopy(cd_input); kt_input["classical_dynamics"]["fewest_switches"] = {}
        lz_input = cp.deepcopy(cd_input); lz_input["classical_dynamics"]["landau_zener"]    = {}

        kt_input["classical_dynamics"]["time_derivative_coupling"] = "kappa"

        fs_tdc_input = cp.deepcopy(fs_input); fs_tdc_input["classical_dynamics"]["trajectories"] = 1;
        kt_tdc_input = cp.deepcopy(kt_input); kt_tdc_input["classical_dynamics"]["trajectories"] = 1;

        fs_tdc_input["classical_dynamics"]["write"]["time_derivative_coupling"] = f"KTSH_COUPLING_TDC_{potential}_FSSH.mat"
        kt_tdc_input["classical_dynamics"]["write"]["time_derivative_coupling"] = f"KTSH_COUPLING_TDC_{potential}_KTSH.mat"

        fs_input["classical_dynamics"]["write"]["population_mean"] = f"KTSH_COUPLING_POPULATION_MEAN_{potential}_FSSH.mat"
        kt_input["classical_dynamics"]["write"]["population_mean"] = f"KTSH_COUPLING_POPULATION_MEAN_{potential}_KTSH.mat"
        lz_input["classical_dynamics"]["write"]["population_mean"] = f"KTSH_COUPLING_POPULATION_MEAN_{potential}_LZSH.mat"

        open(f"input_KTSH_COUPLING_dp_{potential}.json", "w").write(js.dumps(dp_input, indent=4)); os.system(f"acorn input_KTSH_COUPLING_dp_{potential}.json")
        open(f"input_KTSH_COUPLING_ap_{potential}.json", "w").write(js.dumps(ap_input, indent=4)); os.system(f"acorn input_KTSH_COUPLING_ap_{potential}.json")
        open(f"input_KTSH_COUPLING_qd_{potential}.json", "w").write(js.dumps(qd_input, indent=4)); os.system(f"acorn input_KTSH_COUPLING_qd_{potential}.json")
        open(f"input_KTSH_COUPLING_fs_{potential}.json", "w").write(js.dumps(fs_input, indent=4)); os.system(f"acorn input_KTSH_COUPLING_fs_{potential}.json")
        open(f"input_KTSH_COUPLING_kt_{potential}.json", "w").write(js.dumps(kt_input, indent=4)); os.system(f"acorn input_KTSH_COUPLING_kt_{potential}.json")
        open(f"input_KTSH_COUPLING_lz_{potential}.json", "w").write(js.dumps(lz_input, indent=4)); os.system(f"acorn input_KTSH_COUPLING_lz_{potential}.json")

        open(f"input_KTSH_COUPLING_fs_tdc_{potential}.json", "w").write(js.dumps(fs_tdc_input, indent=4)); os.system(f"acorn input_KTSH_COUPLING_fs_tdc_{potential}.json")
        open(f"input_KTSH_COUPLING_kt_tdc_{potential}.json", "w").write(js.dumps(kt_tdc_input, indent=4)); os.system(f"acorn input_KTSH_COUPLING_kt_tdc_{potential}.json")

    os.system(f'lines.py KTSH_COUPLING_POTENTIAL_ADIABATIC_{DS[0]}.mat:0,3 KTSH_COUPLING_TDC_{DS[0]}_FSSH.mat:1 KTSH_COUPLING_TDC_{DS[0]}_KTSH.mat:1 KTSH_COUPLING_POPULATION_{DS[0]}_EXACT.mat:0 KTSH_COUPLING_POPULATION_MEAN_{DS[0]}_FSSH.mat:0 KTSH_COUPLING_POPULATION_MEAN_{DS[0]}_KTSH.mat:0 KTSH_COUPLING_POPULATION_MEAN_{DS[0]}_LZSH.mat:0 KTSH_COUPLING_POTENTIAL_ADIABATIC_{DS[1]}.mat:0,3 KTSH_COUPLING_TDC_{DS[1]}_FSSH.mat:2 KTSH_COUPLING_TDC_{DS[1]}_KTSH.mat:1 KTSH_COUPLING_POPULATION_{DS[1]}_EXACT.mat:0 KTSH_COUPLING_POPULATION_MEAN_{DS[1]}_FSSH.mat:0 KTSH_COUPLING_POPULATION_MEAN_{DS[1]}_KTSH.mat:0 KTSH_COUPLING_POPULATION_MEAN_{DS[1]}_LZSH.mat:0 --dpi 256 --figsize 8 16 --subplots 231 231 232 232 233 233 233 233 234 234 235 235 236 236 236 236 --legends every "S{dollar}_0{dollar}" "S{dollar}_1{dollar}" "T{dollar}_{{0,1}}{dollar} (NPI)" "T{dollar}_{{0,1}}{dollar} ({dollar}{esc}kappa{dollar})" "S{dollar}_0{dollar} (EXACT)" "S{dollar}_0{dollar} (FSSH)" "S{dollar}_0{dollar} ({dollar}{esc}kappa{dollar}TSH)" "S{dollar}_0{dollar} (LZSH)" "S{dollar}_0{dollar}" "S{dollar}_1{dollar}" "T{dollar}_{{0,1}}{dollar} (NPI)" "T{dollar}_{{0,1}}{dollar} ({dollar}{esc}kappa{dollar}TDC)" "S{dollar}_0{dollar} (EXACT)" "S{dollar}_0{dollar} (FSSH)" "S{dollar}_0{dollar} ({dollar}{esc}kappa{dollar}TSH)" "S{dollar}_0{dollar} (LZSH)" --xlabel "Coordinate (a.u.)" "Time (a.u.)" "Time (a.u.)" "Coordinate (a.u.)" "Time (a.u.)" "Time (a.u.)" --xlim nan nan {TDCLIM[DS[0]][0]} {TDCLIM[DS[0]][1]} nan nan nan nan {TDCLIM[DS[1]][0]} {TDCLIM[DS[1]][1]} nan nan --ylabel "Energy (a.u.)" "TDC (a.u.)" "Ground State Population" "Energy (a.u.)" "TDC (a.u.)" "Ground State Population" --output MODEL_COUPLING_DS.png')
    os.system(f'lines.py KTSH_COUPLING_POTENTIAL_ADIABATIC_{TS[0]}.mat:0,4,8 KTSH_COUPLING_TDC_{TS[0]}_FSSH.mat:3,2 KTSH_COUPLING_TDC_{TS[0]}_KTSH.mat:1 KTSH_COUPLING_POPULATION_{TS[0]}_EXACT.mat:0 KTSH_COUPLING_POPULATION_MEAN_{TS[0]}_FSSH.mat:0 KTSH_COUPLING_POPULATION_MEAN_{TS[0]}_KTSH.mat:0 KTSH_COUPLING_POPULATION_MEAN_{TS[0]}_LZSH.mat:0 KTSH_COUPLING_POTENTIAL_ADIABATIC_{TS[1]}.mat:0,4,8 KTSH_COUPLING_TDC_{TS[1]}_FSSH.mat:3,2 KTSH_COUPLING_TDC_{TS[1]}_KTSH.mat:1 KTSH_COUPLING_POPULATION_{TS[1]}_EXACT.mat:0 KTSH_COUPLING_POPULATION_MEAN_{TS[1]}_FSSH.mat:0 KTSH_COUPLING_POPULATION_MEAN_{TS[1]}_KTSH.mat:0 KTSH_COUPLING_POPULATION_MEAN_{TS[1]}_LZSH.mat:0 --dpi 256 --figsize 8 16 --subplots 231 231 231 232 232 232 233 233 233 233 234 234 234 235 235 235 236 236 236 236 --legends every "S{dollar}_0{dollar}" "S{dollar}_1{dollar}" "S{dollar}_2{dollar}" "T{dollar}_{{0,1}}{dollar}=T{dollar}_{{1,2}}{dollar} (NPI)" "T{dollar}_{{0,2}}{dollar} (NPI)" "T{dollar}_{{0,1}}{dollar}=T{dollar}_{{0,2}}{dollar}=T{dollar}_{{1,2}}{dollar} ({dollar}{esc}kappa{dollar}TDC)" "S{dollar}_0{dollar} (EXACT)" "S{dollar}_0{dollar} (FSSH)" "S{dollar}_0{dollar} ({dollar}{esc}kappa{dollar}TSH)" "S{dollar}_0{dollar} (LZSH{dollar}_{esc}mathrm{{maxdiff}}{dollar})" "S{dollar}_0{dollar}" "S{dollar}_1{dollar}" "S{dollar}_2{dollar}" "T{dollar}_{{0,1}}{dollar}=T{dollar}_{{1,2}}{dollar} (NPI)" "T{dollar}_{{0,2}}{dollar} (NPI)" "T{dollar}_{{0,1}}{dollar}=T{dollar}_{{0,2}}{dollar}=T{dollar}_{{1,2}}{dollar} ({dollar}{esc}kappa{dollar}TDC)" "S{dollar}_0{dollar} (EXACT)" "S{dollar}_0{dollar} (FSSH)" "S{dollar}_0{dollar} ({dollar}{esc}kappa{dollar}TSH)" "S{dollar}_0{dollar} (LZSH{dollar}_{esc}mathrm{{maxdiff}}{dollar})" --xlabel "Coordinate (a.u.)" "Time (a.u.)" "Time (a.u.)" "Coordinate (a.u.)" "Time (a.u.)" "Time (a.u.)" --xlim nan nan {TDCLIM[TS[0]][0]} {TDCLIM[TS[0]][1]} nan nan nan nan {TDCLIM[TS[1]][0]} {TDCLIM[TS[1]][1]} nan nan --ylabel "Energy (a.u.)" "TDC (a.u.)" "Ground State Population" "Energy (a.u.)" "TDC (a.u.)" "Ground State Population" --output MODEL_COUPLING_TS.png')

# KTSH TESTING ON MODEL POTENTIALS WITH TRIVIAL CROSSINGS ON HO ========================================================================================================================================

if KTSH_HO:

    HAMILTONIAN = {
        "P1" : {
            "dims" : 1, "states" : 2,
            "matrix" : [
                ["0.0025 * (r1 + 1)^2", "0"                  ],
                ["0",                   "0.0025 * (r1 - 1)^2"],
            ]
        },
        "P2" : {
            "dims" : 1, "states" : 3,
            "matrix" : [
                ["0.0025 * r1^2", "0",                          "0"                   ],
                ["0",             "0.01 + 0.0025 * (r1 - 2)^2", "0"                   ],
                ["0",             "0",                          "0.02 + 0.0025 * r1^2"],
            ]
        },
        "P3" : {
            "dims" : 1, "states" : 3,
            "matrix" : [
                ["0.0025 * r1^2", "0",                          "0"                  ],
                ["0",             "0.05 + 0.0025 * (r1 - 2)^2", "0"                  ],
                ["0",             "0",                          "0.1 + 0.0025 * r1^2"],
            ]
        },
        "P4" : {
            "dims" : 1, "states" : 3,
            "matrix" : [
                ["0.0025 * r1^2",             "0.0002 * exp(-((r1-2)^2))",  "0"                    ],
                ["0.0002 * exp(-((r1-2)^2))", "0.01 + 0.0025 * (r1 - 2)^2", "0.0002 * exp(-(r1^2))"],
                ["0",                         "0.0002 * exp(-(r1^2))",      "0.02 + 0.0025 * r1^2" ],
            ]
        }
    }

    IS = {
        "P1" : [1, 1],
        "P2" : [1, 2],
        "P3" : [1, 2],
        "P4" : [1, 2],
    }

    NS = {
        "P1" : 2,
        "P2" : 3,
        "P3" : 3,
        "P4" : 3,
    }

    TDCLIM = {
        "P1" : [775,  975],
        "P2" : [750,  950],
        "P3" : [400,  600],
        "P4" : [750,  950],
    }

    MODELS = ["P1", "P2", "P3", "P4"]

    for potential in MODELS:

        dp_input = {"model_potential" : {"adiabatic" : False, "limits" : [-10, 10], "output" : f"KTSH_HO_POTENTIAL_DIABATIC_{potential}.mat",  "points" : 8192, "hamiltonian" : HAMILTONIAN[potential]}}
        ap_input = {"model_potential" : {"adiabatic" : True,  "limits" : [-10, 10], "output" : f"KTSH_HO_POTENTIAL_ADIABATIC_{potential}.mat", "points" : 8192, "hamiltonian" : HAMILTONIAN[potential]}}

        qd_input = {
            "quantum_dynamics" : {
                "log_intervals" : {
                    "iteration" : 200
                },

                "grid" : {
                    "limits" : [-24, 24], "points" : NGRID
                },

                "hamiltonian" : HAMILTONIAN[potential],

                "initial_conditions" : {
                    "position" : [-7], "momentum" : [0], "gamma" : 2, "state" : IS[potential][0], "mass" : 2000
                },

                "write" : {
                    "population" : f"KTSH_HO_POPULATION_{potential}_EXACT.mat"
                },

                "adiabatic" : True,
                "iterations" : 800,
                "mode" : [0, 1],
                "time_step" : 10
            }
        }

        cd_input = {
            "classical_dynamics" : {
                "log_intervals" : {
                    "iteration" : 20000, "trajectory" : 1000
                },

                "hamiltonian" : HAMILTONIAN[potential],

                "initial_conditions" : {
                    "position_mean" : qd_input["quantum_dynamics"]["initial_conditions"]["position"], "position_std"  : [0.5],
                    "momentum_mean" : qd_input["quantum_dynamics"]["initial_conditions"]["momentum"], "momentum_std"  : [1.0],
                    "state" : [int(i == IS[potential][1]) for i in range(NS[potential])],
                    "mass" : [2000]
                },

                "write" : {},

                "adiabatic" : qd_input["quantum_dynamics"]["adiabatic"],
                "iterations" : qd_input["quantum_dynamics"]["iterations"] * 10,
                "time_step" : qd_input["quantum_dynamics"]["time_step"] / 10,
                "trajectories" : TRAJS,
                "seed" : 6
            }
        }

        fs_input = cp.deepcopy(cd_input); fs_input["classical_dynamics"]["fewest_switches"] = {}
        lz_input = cp.deepcopy(cd_input); lz_input["classical_dynamics"]["landau_zener"]    = {}
        kt_input = cp.deepcopy(cd_input); kt_input["classical_dynamics"]["fewest_switches"] = {}
        lt_input = cp.deepcopy(cd_input); lt_input["classical_dynamics"]["fewest_switches"] = {}

        kt_input["classical_dynamics"]["time_derivative_coupling"] =  "kappa"
        lt_input["classical_dynamics"]["time_derivative_coupling"] = "lambda"

        fs_tdc_input = cp.deepcopy(fs_input); fs_tdc_input["classical_dynamics"]["trajectories"] = 1;
        kt_tdc_input = cp.deepcopy(kt_input); kt_tdc_input["classical_dynamics"]["trajectories"] = 1;
        lt_tdc_input = cp.deepcopy(lt_input); lt_tdc_input["classical_dynamics"]["trajectories"] = 1;

        fs_tdc_input["classical_dynamics"]["iterations"] *= 10; fs_tdc_input["classical_dynamics"]["time_step"] /= 10
        kt_tdc_input["classical_dynamics"]["iterations"] *= 10; kt_tdc_input["classical_dynamics"]["time_step"] /= 10
        lt_tdc_input["classical_dynamics"]["iterations"] *= 10; lt_tdc_input["classical_dynamics"]["time_step"] /= 10

        fs_tdc_input["classical_dynamics"]["write"]["time_derivative_coupling"] = f"KTSH_HO_TDC_{potential}_FSSH.mat"
        kt_tdc_input["classical_dynamics"]["write"]["time_derivative_coupling"] = f"KTSH_HO_TDC_{potential}_KTSH.mat"
        lt_tdc_input["classical_dynamics"]["write"]["time_derivative_coupling"] = f"KTSH_HO_TDC_{potential}_LTSH.mat"

        fs_input["classical_dynamics"]["write"]["population_mean"] = f"KTSH_HO_POPULATION_MEAN_{potential}_FSSH.mat"
        kt_input["classical_dynamics"]["write"]["population_mean"] = f"KTSH_HO_POPULATION_MEAN_{potential}_KTSH.mat"
        lt_input["classical_dynamics"]["write"]["population_mean"] = f"KTSH_HO_POPULATION_MEAN_{potential}_LTSH.mat"
        lz_input["classical_dynamics"]["write"]["population_mean"] = f"KTSH_HO_POPULATION_MEAN_{potential}_LZSH.mat"

        open(f"input_KTSH_HO_dp_{potential}.json", "w").write(js.dumps(dp_input, indent=4)); os.system(f"acorn input_KTSH_HO_dp_{potential}.json")
        open(f"input_KTSH_HO_ap_{potential}.json", "w").write(js.dumps(ap_input, indent=4)); os.system(f"acorn input_KTSH_HO_ap_{potential}.json")
        open(f"input_KTSH_HO_qd_{potential}.json", "w").write(js.dumps(qd_input, indent=4)); os.system(f"acorn input_KTSH_HO_qd_{potential}.json")
        open(f"input_KTSH_HO_fs_{potential}.json", "w").write(js.dumps(fs_input, indent=4)); os.system(f"acorn input_KTSH_HO_fs_{potential}.json")
        open(f"input_KTSH_HO_kt_{potential}.json", "w").write(js.dumps(kt_input, indent=4)); os.system(f"acorn input_KTSH_HO_kt_{potential}.json")
        open(f"input_KTSH_HO_lt_{potential}.json", "w").write(js.dumps(lt_input, indent=4)); os.system(f"acorn input_KTSH_HO_lt_{potential}.json")
        open(f"input_KTSH_HO_lz_{potential}.json", "w").write(js.dumps(lz_input, indent=4)); os.system(f"acorn input_KTSH_HO_lz_{potential}.json")

        open(f"input_KTSH_HO_fs_tdc_{potential}.json", "w").write(js.dumps(fs_tdc_input, indent=4)); os.system(f"acorn input_KTSH_HO_fs_tdc_{potential}.json")
        open(f"input_KTSH_HO_kt_tdc_{potential}.json", "w").write(js.dumps(kt_tdc_input, indent=4)); os.system(f"acorn input_KTSH_HO_kt_tdc_{potential}.json")
        open(f"input_KTSH_HO_lt_tdc_{potential}.json", "w").write(js.dumps(lt_tdc_input, indent=4)); os.system(f"acorn input_KTSH_HO_lt_tdc_{potential}.json")

    os.system(f'lines.py KTSH_HO_POTENTIAL_ADIABATIC_{MODELS[0]}.mat:0,3 KTSH_HO_TDC_{MODELS[0]}_FSSH.mat:1 KTSH_HO_TDC_{MODELS[0]}_KTSH.mat:1 KTSH_HO_TDC_{MODELS[0]}_LTSH.mat:1 KTSH_HO_POPULATION_{MODELS[0]}_EXACT.mat:1 KTSH_HO_POPULATION_MEAN_{MODELS[0]}_FSSH.mat:1 KTSH_HO_POPULATION_MEAN_{MODELS[0]}_KTSH.mat:1 KTSH_HO_POPULATION_MEAN_{MODELS[0]}_LTSH.mat:1 KTSH_HO_POPULATION_MEAN_{MODELS[0]}_LZSH.mat:1 KTSH_HO_POTENTIAL_ADIABATIC_{MODELS[1]}.mat:0,4,8 KTSH_HO_TDC_{MODELS[1]}_FSSH.mat:5 KTSH_HO_TDC_{MODELS[1]}_KTSH.mat:5 KTSH_HO_TDC_{MODELS[1]}_LTSH.mat:5 KTSH_HO_POPULATION_{MODELS[1]}_EXACT.mat:2 KTSH_HO_POPULATION_MEAN_{MODELS[1]}_FSSH.mat:2 KTSH_HO_POPULATION_MEAN_{MODELS[1]}_KTSH.mat:2 KTSH_HO_POPULATION_MEAN_{MODELS[1]}_LTSH.mat:2 KTSH_HO_POPULATION_MEAN_{MODELS[1]}_LZSH.mat:2 KTSH_HO_POTENTIAL_ADIABATIC_{MODELS[2]}.mat:0,4,8 KTSH_HO_TDC_{MODELS[2]}_FSSH.mat:5 KTSH_HO_TDC_{MODELS[2]}_KTSH.mat:5 KTSH_HO_TDC_{MODELS[2]}_LTSH.mat:5 KTSH_HO_POPULATION_{MODELS[2]}_EXACT.mat:2 KTSH_HO_POPULATION_MEAN_{MODELS[2]}_FSSH.mat:2 KTSH_HO_POPULATION_MEAN_{MODELS[2]}_KTSH.mat:2 KTSH_HO_POPULATION_MEAN_{MODELS[2]}_LTSH.mat:2 KTSH_HO_POPULATION_MEAN_{MODELS[2]}_LZSH.mat:2 KTSH_HO_POTENTIAL_ADIABATIC_{MODELS[3]}.mat:0,4,8 KTSH_HO_TDC_{MODELS[3]}_FSSH.mat:5 KTSH_HO_TDC_{MODELS[3]}_KTSH.mat:5 KTSH_HO_TDC_{MODELS[3]}_LTSH.mat:5 KTSH_HO_POPULATION_{MODELS[3]}_EXACT.mat:2 KTSH_HO_POPULATION_MEAN_{MODELS[3]}_FSSH.mat:2 KTSH_HO_POPULATION_MEAN_{MODELS[3]}_KTSH.mat:2 KTSH_HO_POPULATION_MEAN_{MODELS[3]}_LTSH.mat:2 KTSH_HO_POPULATION_MEAN_{MODELS[3]}_LZSH.mat:2 --dpi 256 --figsize 12 18 --legends every "S{dollar}_0{dollar}" "S{dollar}_1{dollar}" "T{dollar}_{{1,2}}{dollar} (NPI)" "T{dollar}_{{1,2}}{dollar} ({dollar}{esc}kappa{dollar}TDC)" "T{dollar}_{{1,2}}{dollar} ({dollar}{esc}lambda{dollar}TDC)" "S{dollar}_1{dollar} (EXACT)" "S{dollar}_1{dollar} (FSSH)" "S{dollar}_1{dollar} ({dollar}{esc}kappa{dollar}TSH)" "S{dollar}_1{dollar} ({dollar}{esc}lambda{dollar}TSH)" "S{dollar}_1{dollar} (LZSH{dollar}_{esc}mathrm{{maxdiff}}{dollar})" "S{dollar}_0{dollar}" "S{dollar}_1{dollar}" "S{dollar}_2{dollar}" "T{dollar}_{{1,2}}{dollar} (NPI)" "T{dollar}_{{1,2}}{dollar} ({dollar}{esc}kappa{dollar}TDC)" "T{dollar}_{{1,2}}{dollar} ({dollar}{esc}lambda{dollar}TDC)" "S{dollar}_2{dollar} (EXACT)" "S{dollar}_2{dollar} (FSSH)" "S{dollar}_2{dollar} ({dollar}{esc}kappa{dollar}TSH)" "S{dollar}_2{dollar} ({dollar}{esc}lambda{dollar}TSH)" "S{dollar}_2{dollar} (LZSH{dollar}_{esc}mathrm{{maxdiff}}{dollar})" "S{dollar}_0{dollar}" "S{dollar}_1{dollar}" "S{dollar}_2{dollar}" "T{dollar}_{{1,2}}{dollar} (NPI)" "T{dollar}_{{1,2}}{dollar} ({dollar}{esc}kappa{dollar}TDC)" "T{dollar}_{{1,2}}{dollar} ({dollar}{esc}lambda{dollar}TDC)" "S{dollar}_2{dollar} (EXACT)" "S{dollar}_2{dollar} (FSSH)" "S{dollar}_2{dollar} ({dollar}{esc}kappa{dollar}TSH)" "S{dollar}_2{dollar} ({dollar}{esc}lambda{dollar}TSH)" "S{dollar}_2{dollar} (LZSH{dollar}_{esc}mathrm{{maxdiff}}{dollar})" "S{dollar}_0{dollar}" "S{dollar}_1{dollar}" "S{dollar}_2{dollar}" "T{dollar}_{{1,2}}{dollar} (NPI)" "T{dollar}_{{1,2}}{dollar} ({dollar}{esc}kappa{dollar}TDC)" "T{dollar}_{{1,2}}{dollar} ({dollar}{esc}lambda{dollar}TDC)" "S{dollar}_2{dollar} (EXACT)" "S{dollar}_2{dollar} (FSSH)" "S{dollar}_2{dollar} ({dollar}{esc}kappa{dollar}TSH)" "S{dollar}_2{dollar} ({dollar}{esc}lambda{dollar}TSH)" "S{dollar}_2{dollar} (LZSH{dollar}_{esc}mathrm{{maxdiff}}{dollar})" --subplots 431 431 432 432 432 433 433 433 433 433 434 434 434 435 435 435 436 436 436 436 436 437 437 437 438 438 438 439 439 439 439 439 4-3-10 4-3-10 4-3-10 4-3-11 4-3-11 4-3-11 4-3-12 4-3-12 4-3-12 4-3-12 4-3-12 --xlabel "Coordinate (a.u.)" "Time (a.u.)" "Time (a.u.)" "Coordinate (a.u.)" "Time (a.u.)" "Time (a.u.)" "Coordinate (a.u.)" "Time (a.u.)" "Time (a.u.)" "Coordinate (a.u.)" "Time (a.u.)" "Time (a.u.)" --xlim -3 3 {TDCLIM[MODELS[0]][0]} {TDCLIM[MODELS[0]][1]} nan nan -3 3 {TDCLIM[MODELS[1]][0]} {TDCLIM[MODELS[1]][1]} nan nan -8 8 {TDCLIM[MODELS[2]][0]} {TDCLIM[MODELS[2]][1]} nan nan -3 3 {TDCLIM[MODELS[3]][0]} {TDCLIM[MODELS[3]][1]} nan nan --ylabel "Energy (a.u.)" "TDC (a.u.)" "State Population" "Energy (a.u.)" "TDC (a.u.)" "State Population" "Energy (a.u.)" "TDC (a.u.)" "State Population" "Energy (a.u.)" "TDC (a.u.)" "State Population" --ylim -0.002 0.05 -0.002 0.02 nan nan -0.002 0.05 -0.002 0.02 nan nan -0.02 0.25 -0.002 0.02 nan nan -0.002 0.05 -0.005 0.2 nan nan --output MODEL_TRIVIAL_HO_TS.png')

# GENERAL SURFACE HOPPING TESTING ON MODEL POTENTIALS ========================================================================================================================================

if SH_COMPARE:

    HAMILTONIAN = {
        "P1" : {
            "dims" : 1, "states" : 2,
            "matrix" : [
                ["0.01 * tanh(0.6 * r1)", "0.001 * exp(-(r1^2))"  ],
                ["0.001 * exp(-(r1^2))",  "-0.01 * tanh(0.6 * r1)"],
            ]
        }
    }

    IS = {
        "P1" : [1, 1]
    }

    NS = {
        "P1" : 2
    }

    POTENTIALS = ["P1"]

    for potential in POTENTIALS:

        dp_input = {"model_potential" : {"adiabatic" : False, "limits" : [-16, 16], "output" : f"SH_COMPARE_POTENTIAL_DIABATIC_{potential}.mat",  "points" : 8192, "hamiltonian" : HAMILTONIAN[potential]}}
        ap_input = {"model_potential" : {"adiabatic" : True,  "limits" : [-16, 16], "output" : f"SH_COMPARE_POTENTIAL_ADIABATIC_{potential}.mat", "points" : 8192, "hamiltonian" : HAMILTONIAN[potential]}}

        qd_input = {
            "quantum_dynamics" : {
                "log_intervals" : {
                    "iteration" : 50
                },

                "grid" : {
                    "limits" : [-24, 48], "points" : NGRID
                },

                "hamiltonian" : HAMILTONIAN[potential],

                "initial_conditions" : {
                    "position" : [-10], "momentum" : [15], "gamma" : 2, "state" : IS[potential][0], "mass" : 2000
                },

                "write" : {
                    "population" : f"SH_COMPARE_POPULATION_{potential}_EXACT.mat"
                },

                "adiabatic" : True,
                "iterations" : 300,
                "mode" : [0, 1],
                "time_step" : 10
            }
        }

        cd_input = {
            "classical_dynamics" : {
                "log_intervals" : {
                    "iteration" : 1000, "trajectory" : 1000
                },

                "hamiltonian" : HAMILTONIAN[potential],

                "initial_conditions" : {
                    "position_mean" : qd_input["quantum_dynamics"]["initial_conditions"]["position"], "position_std"  : [0.5],
                    "momentum_mean" : qd_input["quantum_dynamics"]["initial_conditions"]["momentum"], "momentum_std"  : [1.0],
                    "state" : [int(i == IS[potential][1]) for i in range(NS[potential])],
                    "mass" : [qd_input["quantum_dynamics"]["initial_conditions"]["mass"]]
                },

                "write" : {},

                "adiabatic" : qd_input["quantum_dynamics"]["adiabatic"],
                "iterations" : qd_input["quantum_dynamics"]["iterations"] * 10,
                "time_step" : qd_input["quantum_dynamics"]["time_step"] / 10,
                "trajectories" : TRAJS
            }
        }

        fs_input = cp.deepcopy(cd_input); fs_input["classical_dynamics"]["fewest_switches"] = {}
        lz_input = cp.deepcopy(cd_input); lz_input["classical_dynamics"]["landau_zener"]    = {}
        ma_input = cp.deepcopy(cd_input); ma_input["classical_dynamics"]["spin_mapping"]    = {}

        fs_input["classical_dynamics"]["write"]["population_mean"] = f"SH_COMPARE_POPULATION_{potential}_FSSH.mat"
        lz_input["classical_dynamics"]["write"]["population_mean"] = f"SH_COMPARE_POPULATION_{potential}_LZSH.mat"
        ma_input["classical_dynamics"]["write"]["population_mean"] = f"SH_COMPARE_POPULATION_{potential}_MASH.mat"

        open(f"input_SH_COMPARE_dp_{potential}.json", "w").write(js.dumps(dp_input, indent=4)); os.system(f"acorn input_SH_COMPARE_dp_{potential}.json")
        open(f"input_SH_COMPARE_ap_{potential}.json", "w").write(js.dumps(ap_input, indent=4)); os.system(f"acorn input_SH_COMPARE_ap_{potential}.json")
        open(f"input_SH_COMPARE_qd_{potential}.json", "w").write(js.dumps(qd_input, indent=4)); os.system(f"acorn input_SH_COMPARE_qd_{potential}.json")
        open(f"input_SH_COMPARE_fs_{potential}.json", "w").write(js.dumps(fs_input, indent=4)); os.system(f"acorn input_SH_COMPARE_fs_{potential}.json")
        open(f"input_SH_COMPARE_lz_{potential}.json", "w").write(js.dumps(lz_input, indent=4)); os.system(f"acorn input_SH_COMPARE_lz_{potential}.json")
        open(f"input_SH_COMPARE_ma_{potential}.json", "w").write(js.dumps(ma_input, indent=4)); os.system(f"acorn input_SH_COMPARE_ma_{potential}.json")

    os.system(f'lines.py SH_COMPARE_POTENTIAL_ADIABATIC_{POTENTIALS[0]}.mat:0,3 SH_COMPARE_POPULATION_{POTENTIALS[0]}_EXACT.mat:0 SH_COMPARE_POPULATION_{POTENTIALS[0]}_FSSH.mat:0 SH_COMPARE_POPULATION_{POTENTIALS[0]}_LZSH.mat:0 SH_COMPARE_POPULATION_{POTENTIALS[0]}_MASH.mat:0 --dpi 256 --figsize 6 16 --subplots 121 121 122 122 122 122 --legends every "S{dollar}_0{dollar}" "S{dollar}_1{dollar}" "S{dollar}_0{dollar} (EXACT)" "S{dollar}_0{dollar} (FSSH)" "S{dollar}_0{dollar} (LZSH)" "S{dollar}_0{dollar} (MASH)" --xlabel "Coordinate (a.u.)" "Time (a.u.)" --ylabel "Energy (a.u.)" "Ground State Population" --output SH_COMPARE.png')

# LZSH FOR THREE STATE CROSSING COMPARE ======================================================================================================================================================

if LZ_COMPARE:

    HAMILTONIAN = {
        "P1" : {
            "dims" : 1, "states" : 3,
            "matrix" : [
                ["sgn(r1) * 0.01 * (1 - exp(-(1.6 * sgn(r1) * r1)))", "0.003 * exp(-(r1^2))", "0"                                                 ],
                ["0.003 * exp(-(r1^2))",                              "0",                    "0.003 * exp(-(r1^2))"                              ],
                ["0",                                                 "0.003 * exp(-(r1^2))", "-sgn(r1) * 0.01 * (1 - exp(-(1.6 * sgn(r1) * r1)))"],
            ]
        },
        "P2" : {
            "dims" : 1, "states" : 3,
            "matrix" : [
                ["sgn(r1) * 0.01 * (1 - exp(-(1.6 * sgn(r1) * r1)))", "0", "0.004 * exp(-(r1^2))"                              ],
                ["0",                                                 "0", "0"                                                 ],
                ["0.004 * exp(-(r1^2))",                              "0", "-sgn(r1) * 0.01 * (1 - exp(-(1.6 * sgn(r1) * r1)))"],
            ]
        },
        "P3" : {
            "dims" : 1, "states" : 3,
            "matrix" : [
                ["0.0025 * r1^2", "0",                          "0"                   ],
                ["0",             "0.01 + 0.0025 * (r1 - 2)^2", "0"                   ],
                ["0",             "0",                          "0.02 + 0.0025 * r1^2"],
            ]
        }
    }

    IS = {
        "P1" : [2, 2],
        "P2" : [2, 2],
        "P3" : [1, 2],
    }

    NS = {
        "P1" : 3,
        "P2" : 3,
        "P3" : 3,
    }

    X = {
        "P1" : -5,
        "P2" : -5,
        "P3" : -2,
    }

    P = {
        "P1" : 15,
        "P2" : 15,
        "P3" :  0,
    }

    POTXLIM = {
        "P1" : [-10, 10],
        "P2" : [-10, 10],
        "P3" : [-3,  3 ],
    }

    POTYLIM = {
        "P1" : [-0.012, 0.012],
        "P2" : [-0.012, 0.012],
        "P3" : [-0.002, 0.050],
    }

    ITERS = {
        "P1" : 150,
        "P2" : 150,
        "P3" : 400,
    }

    MODELS = ["P1", "P2", "P3"]

    for potential in MODELS:

        dp_input = {"model_potential" : {"adiabatic" : False, "limits" : [-16, 16], "output" : f"LZ_COMPARE_POTENTIAL_DIABATIC_{potential}.mat",  "points" : 8192, "hamiltonian" : HAMILTONIAN[potential]}}
        ap_input = {"model_potential" : {"adiabatic" : True,  "limits" : [-16, 16], "output" : f"LZ_COMPARE_POTENTIAL_ADIABATIC_{potential}.mat", "points" : 8192, "hamiltonian" : HAMILTONIAN[potential]}}

        qd_input = {
            "quantum_dynamics" : {
                "log_intervals" : {
                    "iteration" : 50
                },

                "grid" : {
                    "limits" : [-24, 48], "points" : NGRID
                },

                "hamiltonian" : HAMILTONIAN[potential],

                "initial_conditions" : {
                    "position" : [X[potential]], "momentum" : [P[potential]], "gamma" : 2, "state" : IS[potential][0], "mass" : 2000
                },

                "write" : {
                    "population" : f"LZ_COMPARE_POPULATION_{potential}_EXACT.mat"
                },

                "adiabatic" : True,
                "iterations" : ITERS[potential],
                "mode" : [0, 1],
                "time_step" : 10
            }
        }

        cd_input = {
            "classical_dynamics" : {
                "log_intervals" : {
                    "iteration" : 1000, "trajectory" : 1000
                },

                "hamiltonian" : HAMILTONIAN[potential],

                "initial_conditions" : {
                    "position_mean" : qd_input["quantum_dynamics"]["initial_conditions"]["position"], "position_std"  : [0.5],
                    "momentum_mean" : qd_input["quantum_dynamics"]["initial_conditions"]["momentum"], "momentum_std"  : [1.0],
                    "state" : [int(i == IS[potential][1]) for i in range(NS[potential])],
                    "mass" : [qd_input["quantum_dynamics"]["initial_conditions"]["mass"]]
                },

                "write" : {},

                "adiabatic" : qd_input["quantum_dynamics"]["adiabatic"],
                "iterations" : qd_input["quantum_dynamics"]["iterations"] * 10,
                "time_step" : qd_input["quantum_dynamics"]["time_step"] / 10,
                "trajectories" : TRAJS
            }
        }

        lz_prob_input    = cp.deepcopy(cd_input); lz_prob_input   ["classical_dynamics"]["landau_zener"] = {"mode" : None     }
        lz_nearest_input = cp.deepcopy(cd_input); lz_nearest_input["classical_dynamics"]["landau_zener"] = {"mode" : "nearest"}
        lz_maxprob_input = cp.deepcopy(cd_input); lz_maxprob_input["classical_dynamics"]["landau_zener"] = {"mode" : "maxprob"}
        lz_maxdiff_input = cp.deepcopy(cd_input); lz_maxdiff_input["classical_dynamics"]["landau_zener"] = {"mode" : "maxdiff"}

        lz_prob_input   ["classical_dynamics"]["write"]["population_mean"] = f"LZ_COMPARE_POPULATION_{potential}_LZSH-PROB.mat"
        lz_nearest_input["classical_dynamics"]["write"]["population_mean"] = f"LZ_COMPARE_POPULATION_{potential}_LZSH-NEAREST.mat"
        lz_maxprob_input["classical_dynamics"]["write"]["population_mean"] = f"LZ_COMPARE_POPULATION_{potential}_LZSH-MAXPROB.mat"
        lz_maxdiff_input["classical_dynamics"]["write"]["population_mean"] = f"LZ_COMPARE_POPULATION_{potential}_LZSH-MAXDIFF.mat"

        open(f"input_LZ_COMPARE_dp_{potential}.json", "w").write(js.dumps(dp_input, indent=4)); os.system(f"acorn input_LZ_COMPARE_dp_{potential}.json")
        open(f"input_LZ_COMPARE_ap_{potential}.json", "w").write(js.dumps(ap_input, indent=4)); os.system(f"acorn input_LZ_COMPARE_ap_{potential}.json")
        open(f"input_LZ_COMPARE_qd_{potential}.json", "w").write(js.dumps(qd_input, indent=4)); os.system(f"acorn input_LZ_COMPARE_qd_{potential}.json")

        open(f"input_LZ_COMPARE_lz-prob_{potential}.json",    "w").write(js.dumps(lz_prob_input,    indent=4)); os.system(f"acorn input_LZ_COMPARE_lz-prob_{potential}.json"   )
        open(f"input_LZ_COMPARE_lz-nearest_{potential}.json", "w").write(js.dumps(lz_nearest_input, indent=4)); os.system(f"acorn input_LZ_COMPARE_lz-nearest_{potential}.json")
        open(f"input_LZ_COMPARE_lz-maxprob_{potential}.json", "w").write(js.dumps(lz_maxprob_input, indent=4)); os.system(f"acorn input_LZ_COMPARE_lz-maxprob_{potential}.json")
        open(f"input_LZ_COMPARE_lz-maxdiff_{potential}.json", "w").write(js.dumps(lz_maxdiff_input, indent=4)); os.system(f"acorn input_LZ_COMPARE_lz-maxdiff_{potential}.json")

    os.system(f'lines.py LZ_COMPARE_POTENTIAL_ADIABATIC_{MODELS[0]}.mat:0,4,8 LZ_COMPARE_POPULATION_{MODELS[0]}_EXACT.mat LZ_COMPARE_POPULATION_{MODELS[0]}_LZSH-PROB.mat LZ_COMPARE_POPULATION_{MODELS[0]}_LZSH-NEAREST.mat LZ_COMPARE_POPULATION_{MODELS[0]}_LZSH-MAXPROB.mat LZ_COMPARE_POPULATION_{MODELS[0]}_LZSH-MAXDIFF.mat LZ_COMPARE_POTENTIAL_ADIABATIC_{MODELS[1]}.mat:0,4,8 LZ_COMPARE_POPULATION_{MODELS[1]}_EXACT.mat LZ_COMPARE_POPULATION_{MODELS[1]}_LZSH-PROB.mat LZ_COMPARE_POPULATION_{MODELS[1]}_LZSH-NEAREST.mat LZ_COMPARE_POPULATION_{MODELS[1]}_LZSH-MAXPROB.mat LZ_COMPARE_POPULATION_{MODELS[1]}_LZSH-MAXDIFF.mat LZ_COMPARE_POTENTIAL_ADIABATIC_{MODELS[2]}.mat:0,4,8 LZ_COMPARE_POPULATION_{MODELS[2]}_EXACT.mat LZ_COMPARE_POPULATION_{MODELS[2]}_LZSH-PROB.mat LZ_COMPARE_POPULATION_{MODELS[2]}_LZSH-NEAREST.mat LZ_COMPARE_POPULATION_{MODELS[2]}_LZSH-MAXPROB.mat LZ_COMPARE_POPULATION_{MODELS[2]}_LZSH-MAXDIFF.mat --figsize 9 24 --subplots 361 361 361 362 362 362 363 363 363 364 364 364 365 365 365 366 366 366 367 367 367 368 368 368 369 369 369 3-6-10 3-6-10 3-6-10 3-6-11 3-6-11 3-6-11 3-6-12 3-6-12 3-6-12 3-6-13 3-6-13 3-6-13 3-6-14 3-6-14 3-6-14 3-6-15 3-6-15 3-6-15 3-6-16 3-6-16 3-6-16 3-6-17 3-6-17 3-6-17 3-6-18 3-6-18 3-6-18 --title "Potential" "Exact" "Normalized Probability" "Nearest State" "Maximum Probability" "Maximum 2. Derivative" --xlim {POTXLIM[MODELS[0]][0]} {POTXLIM[MODELS[0]][1]} nan nan nan nan nan nan nan nan nan nan {POTXLIM[MODELS[1]][0]} {POTXLIM[MODELS[1]][1]} nan nan nan nan nan nan nan nan nan nan {POTXLIM[MODELS[2]][0]} {POTXLIM[MODELS[2]][1]} nan nan nan nan nan nan nan nan nan nan --ylim {POTYLIM[MODELS[0]][0]} {POTYLIM[MODELS[0]][1]} nan nan nan nan nan nan nan nan nan nan {POTYLIM[MODELS[1]][0]} {POTYLIM[MODELS[1]][1]} nan nan nan nan nan nan nan nan nan nan {POTYLIM[MODELS[2]][0]} {POTYLIM[MODELS[2]][1]} nan nan nan nan nan nan nan nan nan nan --output LZ_COMPARE.png')

# DYNAMICS ON URACIL VC MODEL ================================================================================================================================================================

MAITRA_MCTDH_D0_X = [0.00000000000, 32.40831373777, 87.25315237093, 149.57683263588, 191.95693521605, 234.33703779621, 289.18187642937, 334.05492622013, 353.99850390492, 383.91387043210, 413.82923695927, 456.20933953944, 486.12470606662, 513.54712538319, 530.99775585738, 560.91312238456, 585.84259449054, 613.26501380712, 675.58869407207, 732.92647991582, 785.27837133838, 880.01036534111, 952.30583444845, 999.67183144981, 1074.46024776776, 1134.29098082211, 1181.65697782348, 1243.98065808843, 1348.68444093355, 1438.43054051507, 1530.66958730720, 1607.95095083574, 1677.75347273249, 1752.54188905043, 1832.31619978957, 1897.13282726512, 2001.83661011023, 2061.66734316459, 2141.44165390372, 2193.79354532628, 2261.10312001243, 2325.91974748798, 2398.21521659532, 2477.98952733446]
MAITRA_MCTDH_D0_Y = [0.00911213548,  0.00998573466,  0.01997146932,   0.02710413694,   0.03708987161,   0.05848787446,   0.07417974322,   0.09415121255,   0.11840228245,   0.14835948644,   0.17403708987,   0.20542082738,   0.23395149786,   0.27246790299,   0.29671897289,   0.33380884450,   0.36376604850,   0.38373751783,   0.41084165477,   0.42938659058,   0.46932952924,   0.46504992867,   0.48074179743,   0.50499286733,    0.51925820256,    0.53780313837,    0.56062767475,    0.56776034236,    0.58059914407,    0.58915834522,    0.59486447931,    0.61198288159,    0.61198288159,    0.61198288159,    0.60485021398,    0.60485021398,    0.63052781740,    0.63195435092,    0.64907275320,    0.66333808844,    0.65905848787,    0.68616262482,    0.69900142653,    0.70613409415]
MAITRA_MCTDH_D1_X = [4.127258693094842, 59.845251049875216, 134.13590752558238, 202.23567596164727, 264.1445563580699, 328.11706610103994, 406.534981269842, 497.3346725179285, 614.9615452711315, 744.9701941036191, 879.1061016292014, 1015.3056385013311, 1104.0417004028702, 1192.7777623044094, 1341.3590752558237, 1452.7950599693845, 1543.594751217471, 1615.8217783466307, 1700.430581555075, 1787.1030141100666, 1871.711817318511, 1964.575137913145, 2040.9294237353995, 2117.283709557654, 2216.3379181919304, 2303.010350746922, 2373.1737485295344, 2433.0189995794094, 2482.5461038965477]
MAITRA_MCTDH_D1_Y = [0.048134777376654635, 0.048134777376654635, 0.056558363417569195, 0.07220216606498195, 0.09386281588447654, 0.12154031287605295, 0.1552346570397112, 0.16125150421179302, 0.15042117930204574, 0.1516245487364621, 0.16606498194945848, 0.1624548736462094, 0.15643802647412755, 0.15282791817087846, 0.16726835138387486, 0.1768953068592058, 0.18892900120336945, 0.20336943441636582, 0.21780986762936222, 0.23104693140794225, 0.24187725631768953, 0.24428399518652227, 0.23104693140794225, 0.2226233453670277, 0.21419975932611313, 0.2009626955475331, 0.18772563176895307, 0.18531889290012035, 0.18411552346570398]
MAITRA_MCTDH_D2_X = [0, 35.08169889130616, 96.99057928772879, 163.02671837724628, 202.23567596164727, 253.8264096253328, 297.16262590282867, 334.3079541406822, 359.07150629925127, 377.6441704181781, 410.6622399629368, 443.68030950769554, 466.38023231971715, 489.0801551317388, 517.9709659834027, 546.8617768350666, 567.4980703005408, 594.3252518056573, 633.5342093900583, 668.6159082813645, 722.2702712915974, 796.5609277673045, 889.4242483619386, 980.223939610025, 1056.5782254322796, 1116.423476482155, 1192.7777623044094, 1275.3229361663064, 1345.4863339489186, 1428.0315078108154, 1527.0857164450917, 1599.3127435742513, 1675.667029396506, 1756.1485739118555, 1834.5664890806574, 1923.3025509821964, 2026.4840183095675, 2115.2200802111065, 2218.4015475384776, 2325.7102735589438, 2393.8100419950083, 2490.800621282737]
MAITRA_MCTDH_D2_Y = [0.9386281588447654, 0.9338146811070999, 0.9229843561973526, 0.9097472924187726, 0.8820697954271962, 0.8495788206979543, 0.8074608904933815, 0.7773766546329723, 0.7400722021660651, 0.7027677496991577, 0.6714801444043321, 0.6353790613718412, 0.601684717208183, 0.5716004813477737, 0.5415162454873647, 0.5114320096269555, 0.4873646209386282, 0.46690734055354993, 0.4452466907340554, 0.42478941034897716, 0.40312876052948254, 0.37906137184115524, 0.35499398315282793, 0.3309265944645006, 0.3140794223826715, 0.2900120336943442, 0.2719614921780987, 0.25872442839951865, 0.246690734055355, 0.23345367027677497, 0.21299638989169675, 0.19614921780986763, 0.17087845968712395, 0.16004813477737667, 0.15042117930204574, 0.1407942238267148, 0.12996389891696752, 0.12154031287605295, 0.11793020457280386, 0.11191335740072203, 0.11311672683513839, 0.10830324909747292]
MAITRA_SH_D0_X = [2.2611508427985414, 70.09567612675478, 133.40789972511394, 196.7201233234731, 262.2934977646308, 327.8668722057885, 384.39564327575204, 438.663263502917, 477.10282783049223, 531.3704480576572, 574.3323140708295, 635.3833868263902, 732.6128730667274, 814.0143034074749, 908.9826388050136, 999.4286725169553, 1092.1358570716955, 1175.7984382552415, 1297.9005837663628, 1399.652371692297, 1519.4933663606198, 1623.5063051293528, 1734.3026964264814, 1827.0098809812214, 1935.5451214355514, 2037.2969093614859, 2145.8321498158157, 2263.41199364134, 2371.94723409567, 2478.2213237072015]
MAITRA_SH_D0_Y = [0.006587615283267458, 0.015810276679841896, 0.017127799736495388, 0.019762845849802372, 0.030303030303030304, 0.05138339920948617, 0.08695652173913043, 0.11857707509881424, 0.16600790513833993, 0.20816864295125165, 0.23978919631093545, 0.2621870882740448, 0.28194993412384717, 0.2872200263504611, 0.2951251646903821, 0.30171277997364954, 0.3227931488801054, 0.3412384716732543, 0.3425559947299078, 0.34914361001317523, 0.36363636363636365, 0.38998682476943347, 0.4057971014492754, 0.4189723320158103, 0.43610013175230566, 0.4505928853754941, 0.466403162055336, 0.48484848484848486, 0.5085638998682477, 0.5230566534914362]
MAITRA_SH_D1_X = [2.2611508427985414, 119.8409946683227, 201.2424250090702, 359.52298400496807, 474.8416769876937, 587.8992191276208, 667.0394986255698, 755.2243814947128, 838.8869626782589, 924.8106947046034, 1019.7790301021422, 1107.9639129712853, 1196.1487958404284, 1252.677566910392, 1343.1236006223337, 1406.4358242206927, 1530.7991205746125, 1652.9012660857338, 1768.2199590684593, 1933.283970592753, 2062.1695686322696, 2166.1825074010026, 2247.58393774175, 2313.157312182908, 2356.11917819608, 2480.48247455]
MAITRA_SH_D1_Y = [0.048748353096179184, 0.05270092226613966, 0.0619235836627141, 0.07905138339920949, 0.10408432147562582, 0.11198945981554678, 0.10408432147562582, 0.11067193675889328, 0.11857707509881424, 0.12779973649538867, 0.12911725955204217, 0.12648221343873517, 0.12648221343873517, 0.13570487483530963, 0.16469038208168643, 0.18050065876152832, 0.18181818181818182, 0.19235836627140976, 0.19367588932806326, 0.19367588932806326, 0.1857707509881423, 0.17918313570487485, 0.17786561264822134, 0.1857707509881423, 0.1857707509881423, 0.19499341238471674]
MAITRA_SH_D2_X = [6.783452528395625, 151.49710646750228, 230.63738596545122, 323.34457052019144, 438.663263502917, 544.9373531144485, 608.2495767128077, 694.1733087391522, 748.4409289663172, 832.1035101498633, 918.0272421762078, 1015.2567284165451, 1092.1358570716955, 1189.3653433120328, 1286.5948295523701, 1404.1746733778941, 1505.9264613038285, 1625.7674559721513, 1718.4746405268916, 1824.7487301384228, 1944.5897248067456, 2053.1249652610754, 2150.3544515014128, 2267.934295326937, 2367.424932410073, 2475.960172864403]
MAITRA_SH_D2_Y = [0.9393939393939394, 0.927536231884058, 0.9051383399209486, 0.8695652173913044, 0.7839262187088274, 0.6837944664031621, 0.6403162055335969, 0.6047430830039526, 0.5823451910408433, 0.5678524374176548, 0.5625823451910409, 0.5652173913043479, 0.552042160737813, 0.5256916996047432, 0.5098814229249012, 0.47167325428194995, 0.4545454545454546, 0.4308300395256917, 0.40052700922266143, 0.383399209486166, 0.3715415019762846, 0.3544137022397892, 0.3544137022397892, 0.33201581027667987, 0.3069828722002635, 0.2766798418972332]

if URACIL_LVC:

    omega = {
        8  : np.array([734, 771, 1193, 1383, 1406, 1673, 1761, 1794]) / 8065.54429 / 27.211324570273
    }

    ap_input_8_10 = {"model_potential" : {"adiabatic" : True, "limits" : [-2.0, 2.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q10.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless8D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 0]}}
    ap_input_8_12 = {"model_potential" : {"adiabatic" : True, "limits" : [-2.0, 2.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q12.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless8D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 1]}}
    ap_input_8_18 = {"model_potential" : {"adiabatic" : True, "limits" : [-2.0, 2.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q18.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless8D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 2]}}
    ap_input_8_20 = {"model_potential" : {"adiabatic" : True, "limits" : [-3.5, 3.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q20.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless8D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 3]}}
    ap_input_8_21 = {"model_potential" : {"adiabatic" : True, "limits" : [-4.5, 4.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q21.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless8D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 4]}}
    ap_input_8_24 = {"model_potential" : {"adiabatic" : True, "limits" : [-3.5, 3.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q24.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless8D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 5]}}
    ap_input_8_25 = {"model_potential" : {"adiabatic" : True, "limits" : [-4.5, 4.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q25.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless8D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 6]}}
    ap_input_8_26 = {"model_potential" : {"adiabatic" : True, "limits" : [-4.5, 4.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q26.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless8D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 7]}}

    ap_input_12_3  = {"model_potential" : {"adiabatic" : True, "limits" : [-4.5, 4.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q3.mat",  "points" : 4096, "hamiltonian" : {"name" : "uracilDimless12D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  0]}}
    ap_input_12_7  = {"model_potential" : {"adiabatic" : True, "limits" : [-3.0, 3.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q7.mat",  "points" : 4096, "hamiltonian" : {"name" : "uracilDimless12D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  1]}}
    ap_input_12_10 = {"model_potential" : {"adiabatic" : True, "limits" : [-2.0, 2.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q10.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless12D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  2]}}
    ap_input_12_11 = {"model_potential" : {"adiabatic" : True, "limits" : [-3.0, 3.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q11.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless12D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  3]}}
    ap_input_12_12 = {"model_potential" : {"adiabatic" : True, "limits" : [-2.0, 2.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q12.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless12D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  4]}}
    ap_input_12_18 = {"model_potential" : {"adiabatic" : True, "limits" : [-2.0, 2.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q18.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless12D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  5]}}
    ap_input_12_19 = {"model_potential" : {"adiabatic" : True, "limits" : [-3.5, 3.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q19.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless12D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  6]}}
    ap_input_12_20 = {"model_potential" : {"adiabatic" : True, "limits" : [-3.5, 3.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q20.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless12D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  7]}}
    ap_input_12_21 = {"model_potential" : {"adiabatic" : True, "limits" : [-4.5, 4.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q21.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless12D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  8]}}
    ap_input_12_24 = {"model_potential" : {"adiabatic" : True, "limits" : [-3.5, 3.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q24.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless12D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  9]}}
    ap_input_12_25 = {"model_potential" : {"adiabatic" : True, "limits" : [-4.5, 4.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q25.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless12D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(12) if i != 10]}}
    ap_input_12_26 = {"model_potential" : {"adiabatic" : True, "limits" : [-4.5, 4.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q26.mat", "points" : 4096, "hamiltonian" : {"name" : "uracilDimless12D_1"}, "constant" : [{"index" : i, "value" : 0} for i in range(12) if i != 11]}}

    ap_input_8_10 ["model_potential"]["constant"][3]["value"] = -4
    ap_input_8_12 ["model_potential"]["constant"][3]["value"] = -4
    ap_input_12_10["model_potential"]["constant"][7]["value"] = -4
    ap_input_12_12["model_potential"]["constant"][7]["value"] = -4

    open(f"input_URACIL_LVC_ap_uracil8D_1_10.json", "w").write(js.dumps(ap_input_8_10, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil8D_1_10.json")
    open(f"input_URACIL_LVC_ap_uracil8D_1_12.json", "w").write(js.dumps(ap_input_8_12, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil8D_1_12.json")
    open(f"input_URACIL_LVC_ap_uracil8D_1_18.json", "w").write(js.dumps(ap_input_8_18, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil8D_1_18.json")
    open(f"input_URACIL_LVC_ap_uracil8D_1_20.json", "w").write(js.dumps(ap_input_8_20, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil8D_1_20.json")
    open(f"input_URACIL_LVC_ap_uracil8D_1_21.json", "w").write(js.dumps(ap_input_8_21, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil8D_1_21.json")
    open(f"input_URACIL_LVC_ap_uracil8D_1_24.json", "w").write(js.dumps(ap_input_8_24, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil8D_1_24.json")
    open(f"input_URACIL_LVC_ap_uracil8D_1_25.json", "w").write(js.dumps(ap_input_8_25, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil8D_1_25.json")
    open(f"input_URACIL_LVC_ap_uracil8D_1_26.json", "w").write(js.dumps(ap_input_8_26, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil8D_1_26.json")

    open(f"input_URACIL_LVC_ap_uracil12D_1_3.json",  "w").write(js.dumps(ap_input_12_3,  indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil12D_1_3.json" )
    open(f"input_URACIL_LVC_ap_uracil12D_1_7.json",  "w").write(js.dumps(ap_input_12_7,  indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil12D_1_7.json" )
    open(f"input_URACIL_LVC_ap_uracil12D_1_10.json", "w").write(js.dumps(ap_input_12_10, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil12D_1_10.json")
    open(f"input_URACIL_LVC_ap_uracil12D_1_11.json", "w").write(js.dumps(ap_input_12_11, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil12D_1_11.json")
    open(f"input_URACIL_LVC_ap_uracil12D_1_12.json", "w").write(js.dumps(ap_input_12_12, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil12D_1_12.json")
    open(f"input_URACIL_LVC_ap_uracil12D_1_18.json", "w").write(js.dumps(ap_input_12_18, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil12D_1_18.json")
    open(f"input_URACIL_LVC_ap_uracil12D_1_19.json", "w").write(js.dumps(ap_input_12_19, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil12D_1_19.json")
    open(f"input_URACIL_LVC_ap_uracil12D_1_20.json", "w").write(js.dumps(ap_input_12_20, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil12D_1_20.json")
    open(f"input_URACIL_LVC_ap_uracil12D_1_21.json", "w").write(js.dumps(ap_input_12_21, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil12D_1_21.json")
    open(f"input_URACIL_LVC_ap_uracil12D_1_24.json", "w").write(js.dumps(ap_input_12_24, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil12D_1_24.json")
    open(f"input_URACIL_LVC_ap_uracil12D_1_25.json", "w").write(js.dumps(ap_input_12_25, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil12D_1_25.json")
    open(f"input_URACIL_LVC_ap_uracil12D_1_26.json", "w").write(js.dumps(ap_input_12_26, indent=4)); os.system(f"acorn input_URACIL_LVC_ap_uracil12D_1_26.json")

    for i in [8]:

        cd_input = {
            "classical_dynamics" : {
                "log_intervals" : {
                    "iteration" : 50, "trajectory" : 1000
                },

                "hamiltonian" : {
                    "name" : f"uracil{i}D_1"
                },

                "initial_conditions" : {
                    "position_mean" : [0 for i in range(i)], "position_std" : list(np.sqrt(0.5 / omega[i])),
                    "momentum_mean" : [0 for i in range(i)], "momentum_std" : list(np.sqrt(omega[i] / 0.5)),
                    "state" : [0.01, 0.05, 0.94, 0], "mass" : [1 for i in range(i)]
                },

                "write" : {},

                "adiabatic" : True,
                "iterations" : 250,
                "time_step" : 10,
                "trajectories" : TRAJS
            }
        }

        fs_input = cp.deepcopy(cd_input); fs_input["classical_dynamics"]["fewest_switches"] = {}
        kt_input = cp.deepcopy(cd_input); kt_input["classical_dynamics"]["fewest_switches"] = {}

        lz_nearest_input = cp.deepcopy(cd_input); lz_nearest_input["classical_dynamics"]["landau_zener"] = {"mode" : "nearest"}
        lz_maxdiff_input = cp.deepcopy(cd_input); lz_maxdiff_input["classical_dynamics"]["landau_zener"] = {"mode" : "maxdiff"}

        kt_input["classical_dynamics"]["time_derivative_coupling"] = "kappa"

        fs_input["classical_dynamics"]["write"]["population_mean"] = f"URACIL_LVC_POPULATION_uracil{i}D_1_FSSH.mat"
        kt_input["classical_dynamics"]["write"]["population_mean"] = f"URACIL_LVC_POPULATION_uracil{i}D_1_KTSH.mat"

        lz_nearest_input["classical_dynamics"]["write"]["population_mean"] = f"URACIL_LVC_POPULATION_uracil{i}D_1_LZSH-NEAREST.mat"
        lz_maxdiff_input["classical_dynamics"]["write"]["population_mean"] = f"URACIL_LVC_POPULATION_uracil{i}D_1_LZSH-MAXDIFF.mat"

        open(f"input_URACIL_LVC_fs_uracil{i}D_1.json", "w").write(js.dumps(fs_input, indent=4)); os.system(f"acorn input_URACIL_LVC_fs_uracil{i}D_1.json")
        open(f"input_URACIL_LVC_kt_uracil{i}D_1.json", "w").write(js.dumps(kt_input, indent=4)); os.system(f"acorn input_URACIL_LVC_kt_uracil{i}D_1.json")

        open(f"input_URACIL_LVC_lz-nearest_uracil{i}D_1.json", "w").write(js.dumps(lz_nearest_input, indent=4)); os.system(f"acorn input_URACIL_LVC_lz-nearest_uracil{i}D_1.json")
        open(f"input_URACIL_LVC_lz-maxdiff_uracil{i}D_1.json", "w").write(js.dumps(lz_maxdiff_input, indent=4)); os.system(f"acorn input_URACIL_LVC_lz-maxdiff_uracil{i}D_1.json")

        np.savetxt("URACIL_LVC_MAITRA_MCTDH_D0.mat", np.array([MAITRA_MCTDH_D0_X, MAITRA_MCTDH_D0_Y]).T, comments="", fmt="%20.14f", header=f"{len(MAITRA_MCTDH_D0_X)} 2")
        np.savetxt("URACIL_LVC_MAITRA_MCTDH_D1.mat", np.array([MAITRA_MCTDH_D1_X, MAITRA_MCTDH_D1_Y]).T, comments="", fmt="%20.14f", header=f"{len(MAITRA_MCTDH_D1_X)} 2")
        np.savetxt("URACIL_LVC_MAITRA_MCTDH_D2.mat", np.array([MAITRA_MCTDH_D2_X, MAITRA_MCTDH_D2_Y]).T, comments="", fmt="%20.14f", header=f"{len(MAITRA_MCTDH_D2_X)} 2")

        os.system(f'lines.py URACIL_LVC_MAITRA_MCTDH_D0.mat URACIL_LVC_POPULATION_uracil{i}D_1_FSSH.mat:0 URACIL_LVC_POPULATION_uracil{i}D_1_KTSH.mat:0 URACIL_LVC_POPULATION_uracil{i}D_1_LZSH-NEAREST.mat:0 URACIL_LVC_POPULATION_uracil{i}D_1_LZSH-MAXDIFF.mat:0 URACIL_LVC_MAITRA_MCTDH_D1.mat URACIL_LVC_POPULATION_uracil{i}D_1_FSSH.mat:1 URACIL_LVC_POPULATION_uracil{i}D_1_KTSH.mat:1 URACIL_LVC_POPULATION_uracil{i}D_1_LZSH-NEAREST.mat:1 URACIL_LVC_POPULATION_uracil{i}D_1_LZSH-MAXDIFF.mat:1 URACIL_LVC_MAITRA_MCTDH_D2.mat URACIL_LVC_POPULATION_uracil{i}D_1_FSSH.mat:2 URACIL_LVC_POPULATION_uracil{i}D_1_KTSH.mat:2 URACIL_LVC_POPULATION_uracil{i}D_1_LZSH-NEAREST.mat:2 URACIL_LVC_POPULATION_uracil{i}D_1_LZSH-MAXDIFF.mat:2 URACIL_LVC_POPULATION_uracil{i}D_1_FSSH.mat:3 URACIL_LVC_POPULATION_uracil{i}D_1_KTSH.mat:3 URACIL_LVC_POPULATION_uracil{i}D_1_LZSH-NEAREST.mat:3 URACIL_LVC_POPULATION_uracil{i}D_1_LZSH-MAXDIFF.mat:3 --figsize 6 28 --legends every "D{dollar}_0{dollar} (MCTDH)" "D{dollar}_0{dollar} (FSSH)" "D{dollar}_0{dollar} ({dollar}{esc}kappa{dollar}TSH)" "D{dollar}_0{dollar} (LZSH{dollar}_{esc}mathrm{{nearest}}{dollar})" "D{dollar}_0{dollar} (LZSH{dollar}_{esc}mathrm{{maxdiff}}{dollar})" "D{dollar}_1{dollar} (MCTDH)" "D{dollar}_1{dollar} (FSSH)" "D{dollar}_1{dollar} ({dollar}{esc}kappa{dollar}TSH)" "D{dollar}_1{dollar} (LZSH{dollar}_{esc}mathrm{{nearest}}{dollar})" "D{dollar}_1{dollar} (LZSH{dollar}_{esc}mathrm{{maxdiff}}{dollar})" "D{dollar}_2{dollar} (MCTDH)" "D{dollar}_2{dollar} (FSSH)" "D{dollar}_2{dollar} ({dollar}{esc}kappa{dollar}TSH)" "D{dollar}_2{dollar} (LZSH{dollar}_{esc}mathrm{{nearest}}{dollar})" "D{dollar}_2{dollar} (LZSH{dollar}_{esc}mathrm{{maxdiff}}{dollar})" "D{dollar}_3{dollar} (FSSH)" "D{dollar}_3{dollar} ({dollar}{esc}kappa{dollar}TSH)" "D{dollar}_3{dollar} (LZSH{dollar}_{esc}mathrm{{nearest}}{dollar})" "D{dollar}_3{dollar} (LZSH{dollar}_{esc}mathrm{{maxdiff}}{dollar})" --xlabel "Time (a.u.)" "Time (a.u.)" "Time (a.u.)" "Time (a.u.)" --ylabel "Population" "Population" "Population" "Population" --subplots 141 141 141 141 141 142 142 142 142 142 143 143 143 143 143 144 144 144 144 --ylim -0.002 1.002 -0.002 1.002 -0.002 1.002 -0.002 1.002 --output POPULATION_uracil{i}D_1.png')

    os.system(f'lines.py URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q18.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q20.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q21.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q24.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q25.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q26.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q10.mat:0,5,10 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q12.mat:0,5,10 --figsize 6 8 --subplots 241 241 241 241 242 242 242 242 243 243 243 243 244 244 244 244 245 245 245 245 246 246 246 246 247 247 247 248 248 248 --xlabel "Q{dollar}_{{18}}{dollar}" "Q{dollar}_{{20}}{dollar}" "Q{dollar}_{{21}}{dollar}" "Q{dollar}_{{24}}{dollar}" "Q{dollar}_{{25}}{dollar}" "Q{dollar}_{{26}}{dollar}" "Q{dollar}_{{10}}{dollar}" "Q{dollar}_{{12}}{dollar}" --output POTENTIAL_ADIABATIC_uracilDimless8D_1.png')
    os.system(f'lines.py URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q3.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q7.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q11.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q18.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q19.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q20.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q21.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q24.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q25.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q26.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q10.mat:0,5,10 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q12.mat:0,5,10 --figsize 9 8 --subplots 341 341 341 341 342 342 342 342 343 343 343 343 344 344 344 344 345 345 345 345 346 346 346 346 347 347 347 347 348 348 348 348 349 349 349 349 3-4-10 3-4-10 3-4-10 3-4-10 3-4-11 3-4-11 3-4-11 3-4-12 3-4-12 3-4-12 --xlabel "Q{dollar}_{{3}}{dollar}" "Q{dollar}_{{7}}{dollar}" "Q{dollar}_{{11}}{dollar}" "Q{dollar}_{{18}}{dollar}" "Q{dollar}_{{19}}{dollar}" "Q{dollar}_{{20}}{dollar}" "Q{dollar}_{{21}}{dollar}" "Q{dollar}_{{24}}{dollar}" "Q{dollar}_{{25}}{dollar}" "Q{dollar}_{{26}}{dollar}" "Q{dollar}_{{10}}{dollar}" "Q{dollar}_{{12}}{dollar}" --output POTENTIAL_ADIABATIC_uracilDimless12D_1.png')
# CLEANUP ====================================================================================================================================================================================

if os.path.exists("research"): sh.rmtree("research")

os.mkdir("research"); os.mkdir("research/source"); os.mkdir("research/output")

for file in gl.glob("*.json"): sh.move(file, "research/source")
for file in gl.glob("*.mat" ): sh.move(file, "research/source")
for file in gl.glob("*.png" ): sh.move(file, "research/output")
