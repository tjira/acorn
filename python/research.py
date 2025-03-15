#!/usr/bin/env python

import copy as cp, glob as gl, shutil as sh, json as js, numpy as np, os

SH_TESTING = 1; URACIL_LVC = 1

# GENERAL SURFACE HOPPING TESTING ON MODEL POTENTIALS ========================================================================================================================================

if SH_TESTING:

    IS = {
        "doubleState1D_2" : [1, 1],
        "tripleState1D_2" : [2, 2]
    }

    NS = {
        "doubleState1D_2" : 2,
        "tripleState1D_2" : 3
    }

    for potential in ["doubleState1D_2", "tripleState1D_2"]:

        dp_input = {"model_potential" : {"adiabatic" : False, "limits" : [-16, 16], "output" : f"POTENTIAL_DIABATIC_{potential}.mat",  "points" : 8192, "potential" : potential}}
        ap_input = {"model_potential" : {"adiabatic" : True,  "limits" : [-16, 16], "output" : f"POTENTIAL_ADIABATIC_{potential}.mat", "points" : 8192, "potential" : potential}}

        qd_input = {
            "quantum_dynamics" : {
                "log_intervals" : {
                    "iteration" : 50
                },

                "grid" : {
                    "limits" : [-24, 48], "points" : 8192
                },

                "initial_conditions" : {
                    "position" : [-10], "momentum" : [15], "gamma" : 2, "state" : IS[potential][0], "mass" : 2000
                },

                "write" : {
                    "population" : f"POPULATION_{potential}_EXACT.mat"
                },

                "adiabatic" : True,
                "iterations" : 300,
                "mode" : [0, 1],
                "potential" : potential,
                "time_step" : 10
            }
        }

        cd_input = {
            "classical_dynamics" : {
                "log_intervals" : {
                    "iteration" : 500, "trajectory" : 100
                },

                "initial_conditions" : {
                    "position_mean" : qd_input["quantum_dynamics"]["initial_conditions"]["position"], "position_std"  : [0.5],
                    "momentum_mean" : qd_input["quantum_dynamics"]["initial_conditions"]["momentum"], "momentum_std"  : [1.0],
                    "state" : [int(i == IS[potential][1]) for i in range(NS[potential])],
                    "mass" : [qd_input["quantum_dynamics"]["initial_conditions"]["mass"]]
                },

                "write" : {},

                "adiabatic" : qd_input["quantum_dynamics"]["adiabatic"],
                "iterations" : qd_input["quantum_dynamics"]["iterations"] * 10,
                "potential" : potential,
                "time_step" : qd_input["quantum_dynamics"]["time_step"] / 10,
                "trajectories" : 10000
            }
        }

        fs_input = cp.deepcopy(cd_input); fs_input["classical_dynamics"]["fewest_switches"] = {}
        kt_input = cp.deepcopy(cd_input); kt_input["classical_dynamics"]["fewest_switches"] = {}
        lz_input = cp.deepcopy(cd_input); lz_input["classical_dynamics"]["landau_zener"]    = {}

        kt_input["classical_dynamics"]["time_derivative_coupling"] = "baeckan"

        fs_input["classical_dynamics"]["write"]["population_mean"] = f"POPULATION_{potential}_FSSH.mat"
        kt_input["classical_dynamics"]["write"]["population_mean"] = f"POPULATION_{potential}_KTSH.mat"
        lz_input["classical_dynamics"]["write"]["population_mean"] = f"POPULATION_{potential}_LZSH.mat"

        open(f"input_dp_{potential}.json", "w").write(js.dumps(dp_input, indent=4)); os.system(f"acorn input_dp_{potential}.json")
        open(f"input_ap_{potential}.json", "w").write(js.dumps(ap_input, indent=4)); os.system(f"acorn input_ap_{potential}.json")
        open(f"input_qd_{potential}.json", "w").write(js.dumps(qd_input, indent=4)); os.system(f"acorn input_qd_{potential}.json")
        open(f"input_fs_{potential}.json", "w").write(js.dumps(fs_input, indent=4)); os.system(f"acorn input_fs_{potential}.json")
        open(f"input_kt_{potential}.json", "w").write(js.dumps(kt_input, indent=4)); os.system(f"acorn input_kt_{potential}.json")
        open(f"input_lz_{potential}.json", "w").write(js.dumps(lz_input, indent=4)); os.system(f"acorn input_lz_{potential}.json")

        os.system(f'lines.py POPULATION_{potential}_EXACT.mat:0 POPULATION_{potential}_FSSH.mat:0 POPULATION_{potential}_KTSH.mat:0 POPULATION_{potential}_LZSH.mat:0 --legend "S\$_0\$ (EXACT)" "S\$_0\$ (FSSH)" "S\$_0\$ (KTSH)" "S\$_0\$ (LZSH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_{potential}.png')

# DYNAMICS ON URACIL VC MODEL ================================================================================================================================================================

if URACIL_LVC:

    omega = {
        8 : np.array([734, 771, 1193, 1383, 1406, 1673, 1761, 1794]) / 8065.54429 / 27.211324570273
    }

    for i in [8]:

        cd_input = {
            "classical_dynamics" : {
                "log_intervals" : {
                    "iteration" : 500, "trajectory" : 100
                },

                "initial_conditions" : {
                    "position_mean" : [0 for i in range(i)], "position_std" : list(np.sqrt(0.5 / omega[i])),
                    "momentum_mean" : [0 for i in range(i)], "momentum_std" : list(np.sqrt(omega[i] / 0.5)),
                    "state" : [0.01, 0.05, 0.94, 0], "mass" : [1 for i in range(i)]
                },

                "write" : {},

                "adiabatic" : True,
                "iterations" : 250,
                "potential" : f"uracil{i}D_1",
                "time_step" : 10,
                "trajectories" : 10000
            }
        }

        fs_input = cp.deepcopy(cd_input); fs_input["classical_dynamics"]["fewest_switches"] = {}
        kt_input = cp.deepcopy(cd_input); kt_input["classical_dynamics"]["fewest_switches"] = {}
        lz_input = cp.deepcopy(cd_input); lz_input["classical_dynamics"]["landau_zener"]    = {}

        kt_input["classical_dynamics"]["time_derivative_coupling"] = "baeckan"

        fs_input["classical_dynamics"]["write"]["population_mean"] = f"POPULATION_uracil{i}D_1_FSSH.mat"
        kt_input["classical_dynamics"]["write"]["population_mean"] = f"POPULATION_uracil{i}D_1_KTSH.mat"
        lz_input["classical_dynamics"]["write"]["population_mean"] = f"POPULATION_uracil{i}D_1_LZSH.mat"

        open(f"input_fs_uracil{i}D_1.json", "w").write(js.dumps(fs_input, indent=4)); os.system(f"acorn input_fs_uracil{i}D_1.json")
        open(f"input_kt_uracil{i}D_1.json", "w").write(js.dumps(kt_input, indent=4)); os.system(f"acorn input_kt_uracil{i}D_1.json")
        open(f"input_lz_uracil{i}D_1.json", "w").write(js.dumps(lz_input, indent=4)); os.system(f"acorn input_lz_uracil{i}D_1.json")

        os.system(f'lines.py external/MAITRA_URACIL_8D_LVC_D0_MCTDH.dat POPULATION_uracil{i}D_1_FSSH.mat:0 POPULATION_uracil{i}D_1_KTSH.mat:0 POPULATION_uracil{i}D_1_LZSH.mat:0 --legend "D\$_0\$ (MCTDH)" "D\$_0\$ (FSSH)" "D\$_0\$ (KTSH)" "D\$_0\$ (LZSH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_uracil{i}D_1.png')

# CLEANUP ====================================================================================================================================================================================

if os.path.exists("research"): sh.rmtree("research")

os.mkdir("research"); os.mkdir("research/source"); os.mkdir("research/output")

for file in gl.glob("*.json"): sh.move(file, "research/source")
for file in gl.glob("*.mat" ): sh.move(file, "research/source")
for file in gl.glob("*.png" ): sh.move(file, "research/output")
