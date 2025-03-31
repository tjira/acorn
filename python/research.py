#!/usr/bin/env python

import copy as cp, glob as gl, shutil as sh, json as js, numpy as np, os

dollar = "$"  if os.name == "nt" else  "\\$"
esc    = "\\" if os.name == "nt" else "\\\\"

KTSH_TESTING = 1; SH_COMPARE = 1; URACIL_LVC = 1

TRAJS = 10000

# KTSH TESTING ON MODEL POTENTIALS ===========================================================================================================================================================

if KTSH_TESTING:

    IS = {
        "tully1D_1"               : [1, 1],
        "tully1D_1_lorentz"       : [1, 1],
        "tripleState1D_4"         : [2, 2],
        "tripleState1D_4_lorentz" : [2, 2]
    }

    NS = {
        "tully1D_1"               : 2,
        "tully1D_1_lorentz"       : 2,
        "tripleState1D_4"         : 3,
        "tripleState1D_4_lorentz" : 3
    }

    TDCLIM = {
        "tully1D_1"               : [1117, 1717],
        "tully1D_1_lorentz"       : [1117, 1717],
        "tripleState1D_4"         : [1117, 1717],
        "tripleState1D_4_lorentz" : [1117, 1717]
    }

    DS, TS = ["tully1D_1", "tully1D_1_lorentz"], ["tripleState1D_4", "tripleState1D_4_lorentz"]

    for potential in DS + TS:

        dp_input = {"model_potential" : {"adiabatic" : False, "limits" : [-16, 16], "output" : f"KTSH_TESTING_POTENTIAL_DIABATIC_{potential}.mat",  "points" : 8192, "potential" : potential}}
        ap_input = {"model_potential" : {"adiabatic" : True,  "limits" : [-16, 16], "output" : f"KTSH_TESTING_POTENTIAL_ADIABATIC_{potential}.mat", "points" : 8192, "potential" : potential}}

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
                    "population" : f"KTSH_TESTING_POPULATION_{potential}_EXACT.mat"
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
                    "iteration" : 1000, "trajectory" : 1000
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
                "trajectories" : TRAJS,
                "seed" : 30
            }
        }

        fs_input = cp.deepcopy(cd_input); fs_input["classical_dynamics"]["fewest_switches"] = {}
        kt_input = cp.deepcopy(cd_input); kt_input["classical_dynamics"]["fewest_switches"] = {}
        lz_input = cp.deepcopy(cd_input); lz_input["classical_dynamics"]["landau_zener"]    = {}

        kt_input["classical_dynamics"]["time_derivative_coupling"] = "baeckan"

        fs_input["classical_dynamics"]["write"]["population_mean"] = f"KTSH_TESTING_POPULATION_MEAN_{potential}_FSSH.mat"
        kt_input["classical_dynamics"]["write"]["population_mean"] = f"KTSH_TESTING_POPULATION_MEAN_{potential}_KTSH.mat"
        lz_input["classical_dynamics"]["write"]["population_mean"] = f"KTSH_TESTING_POPULATION_MEAN_{potential}_LZSH.mat"

        fs_input["classical_dynamics"]["write"]["time_derivative_coupling"] = f"KTSH_TESTING_TDC_{potential}_FSSH.mat"
        kt_input["classical_dynamics"]["write"]["time_derivative_coupling"] = f"KTSH_TESTING_TDC_{potential}_KTSH.mat"

        open(f"input_KTSH_TESTING_dp_{potential}.json", "w").write(js.dumps(dp_input, indent=4)); os.system(f"acorn input_KTSH_TESTING_dp_{potential}.json")
        open(f"input_KTSH_TESTING_ap_{potential}.json", "w").write(js.dumps(ap_input, indent=4)); os.system(f"acorn input_KTSH_TESTING_ap_{potential}.json")
        open(f"input_KTSH_TESTING_qd_{potential}.json", "w").write(js.dumps(qd_input, indent=4)); os.system(f"acorn input_KTSH_TESTING_qd_{potential}.json")
        open(f"input_KTSH_TESTING_fs_{potential}.json", "w").write(js.dumps(fs_input, indent=4)); os.system(f"acorn input_KTSH_TESTING_fs_{potential}.json")
        open(f"input_KTSH_TESTING_kt_{potential}.json", "w").write(js.dumps(kt_input, indent=4)); os.system(f"acorn input_KTSH_TESTING_kt_{potential}.json")
        open(f"input_KTSH_TESTING_lz_{potential}.json", "w").write(js.dumps(lz_input, indent=4)); os.system(f"acorn input_KTSH_TESTING_lz_{potential}.json")

    os.system(f'python python/lines.py KTSH_TESTING_POTENTIAL_ADIABATIC_{DS[0]}.mat:0,3 KTSH_TESTING_TDC_{DS[0]}_FSSH.mat:1 KTSH_TESTING_TDC_{DS[0]}_KTSH.mat:1 KTSH_TESTING_POPULATION_{DS[0]}_EXACT.mat:0 KTSH_TESTING_POPULATION_MEAN_{DS[0]}_FSSH.mat:0 KTSH_TESTING_POPULATION_MEAN_{DS[0]}_KTSH.mat:0 KTSH_TESTING_POPULATION_MEAN_{DS[0]}_LZSH.mat:0 KTSH_TESTING_POTENTIAL_ADIABATIC_{DS[1]}.mat:0,3 KTSH_TESTING_TDC_{DS[1]}_FSSH.mat:1 KTSH_TESTING_TDC_{DS[1]}_KTSH.mat:1 KTSH_TESTING_POPULATION_{DS[1]}_EXACT.mat:0 KTSH_TESTING_POPULATION_MEAN_{DS[1]}_FSSH.mat:0 KTSH_TESTING_POPULATION_MEAN_{DS[1]}_KTSH.mat:0 KTSH_TESTING_POPULATION_MEAN_{DS[1]}_LZSH.mat:0 --dpi 256 --figsize 8 16 --subplots 231 231 232 232 233 233 233 233 234 234 235 235 236 236 236 236 --legends every "S{dollar}_0{dollar}" "S{dollar}_1{dollar}" "T{dollar}_{{0,1}}{dollar} (HST)" "T{dollar}_{{0,1}}{dollar} (BA)" "S{dollar}_0{dollar} (EXACT)" "S{dollar}_0{dollar} (FSSH)" "S{dollar}_0{dollar} ({dollar}{esc}kappa{dollar}TSH)" "S{dollar}_0{dollar} (LZSH)" "S{dollar}_0{dollar}" "S{dollar}_1{dollar}" "T{dollar}_{{0,1}}{dollar} (HST)" "T{dollar}_{{0,1}}{dollar} (BA)" "S{dollar}_0{dollar} (EXACT)" "S{dollar}_0{dollar} (FSSH)" "S{dollar}_0{dollar} ({dollar}{esc}kappa{dollar}TSH)" "S{dollar}_0{dollar} (LZSH)" --xlabel "Coordinate (a.u.)" "Time (a.u.)" "Time (a.u.)" "Coordinate (a.u.)" "Time (a.u.)" "Time (a.u.)" --xlim nan nan {TDCLIM[DS[0]][0]} {TDCLIM[DS[0]][1]} nan nan nan nan {TDCLIM[DS[1]][0]} {TDCLIM[DS[1]][1]} nan nan --ylabel "Energy (a.u.)" "TDC (??)" "Ground State Population" "Energy (a.u.)" "TDC (??)" "Ground State Population" --output MODEL_DS.png')
    os.system(f'python python/lines.py KTSH_TESTING_POTENTIAL_ADIABATIC_{TS[0]}.mat:0,4,8 KTSH_TESTING_TDC_{TS[0]}_FSSH.mat:1,2 KTSH_TESTING_TDC_{TS[0]}_KTSH.mat:1 KTSH_TESTING_POPULATION_{TS[0]}_EXACT.mat:0 KTSH_TESTING_POPULATION_MEAN_{TS[0]}_FSSH.mat:0 KTSH_TESTING_POPULATION_MEAN_{TS[0]}_KTSH.mat:0 KTSH_TESTING_POPULATION_MEAN_{TS[0]}_LZSH.mat:0 KTSH_TESTING_POTENTIAL_ADIABATIC_{TS[1]}.mat:0,4,8 KTSH_TESTING_TDC_{TS[1]}_FSSH.mat:1,2 KTSH_TESTING_TDC_{TS[1]}_KTSH.mat:1 KTSH_TESTING_POPULATION_{TS[1]}_EXACT.mat:0 KTSH_TESTING_POPULATION_MEAN_{TS[1]}_FSSH.mat:0 KTSH_TESTING_POPULATION_MEAN_{TS[1]}_KTSH.mat:0 KTSH_TESTING_POPULATION_MEAN_{TS[1]}_LZSH.mat:0 --dpi 256 --figsize 8 16 --subplots 231 231 231 232 232 232 233 233 233 233 234 234 234 235 235 235 236 236 236 236 --legends every "S{dollar}_0{dollar}" "S{dollar}_1{dollar}" "S{dollar}_2{dollar}" "T{dollar}_{{0,1}}{dollar}=T{dollar}_{{1,2}}{dollar} (HST)" "T{dollar}_{{0,2}}{dollar} (HST)" "T{dollar}_{{0,1}}{dollar}=T{dollar}_{{0,2}}{dollar}=T{dollar}_{{1,2}}{dollar} (BA)" "S{dollar}_0{dollar} (EXACT)" "S{dollar}_0{dollar} (FSSH)" "S{dollar}_0{dollar} ({dollar}{esc}kappa{dollar}TSH)" "S{dollar}_0{dollar} (LZSH)" "S{dollar}_0{dollar}" "S{dollar}_1{dollar}" "S{dollar}_2{dollar}" "T{dollar}_{{0,1}}{dollar}=T{dollar}_{{1,2}}{dollar} (HST)" "T{dollar}_{{0,1}}{dollar} (BA)" "T{dollar}_{{0,1}}{dollar}=T{dollar}_{{0,2}}{dollar}=T{dollar}_{{1,2}}{dollar} (BA)" "S{dollar}_0{dollar} (EXACT)" "S{dollar}_0{dollar} (FSSH)" "S{dollar}_0{dollar} ({dollar}{esc}kappa{dollar}TSH)" "S{dollar}_0{dollar} (LZSH)" --xlabel "Coordinate (a.u.)" "Time (a.u.)" "Time (a.u.)" "Coordinate (a.u.)" "Time (a.u.)" "Time (a.u.)" --xlim nan nan {TDCLIM[TS[0]][0]} {TDCLIM[TS[0]][1]} nan nan nan nan {TDCLIM[TS[1]][0]} {TDCLIM[TS[1]][1]} nan nan --ylabel "Energy (a.u.)" "TDC (??)" "Ground State Population" "Energy (a.u.)" "TDC (??)" "Ground State Population" --output MODEL_TS.png')

# GENERAL SURFACE HOPPING TESTING ON MODEL POTENTIALS ========================================================================================================================================

if SH_COMPARE:

    IS = {
        "doubleState1D_2" : [1, 1]
    }

    NS = {
        "doubleState1D_2" : 2
    }

    POTENTIALS = ["doubleState1D_2"]

    for potential in POTENTIALS:

        dp_input = {"model_potential" : {"adiabatic" : False, "limits" : [-16, 16], "output" : f"SH_COMPARE_POTENTIAL_DIABATIC_{potential}.mat",  "points" : 8192, "potential" : potential}}
        ap_input = {"model_potential" : {"adiabatic" : True,  "limits" : [-16, 16], "output" : f"SH_COMPARE_POTENTIAL_ADIABATIC_{potential}.mat", "points" : 8192, "potential" : potential}}

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
                    "population" : f"SH_COMPARE_POPULATION_{potential}_EXACT.mat"
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
                    "iteration" : 1000, "trajectory" : 1000
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

    os.system(f'python python/lines.py SH_COMPARE_POTENTIAL_ADIABATIC_{POTENTIALS[0]}.mat:0,3 SH_COMPARE_POPULATION_{POTENTIALS[0]}_EXACT.mat:0 SH_COMPARE_POPULATION_{POTENTIALS[0]}_FSSH.mat:0 SH_COMPARE_POPULATION_{POTENTIALS[0]}_LZSH.mat:0 SH_COMPARE_POPULATION_{POTENTIALS[0]}_MASH.mat:0 --dpi 256 --figsize 6 16 --subplots 121 121 122 122 122 122 --legends every "S{dollar}_0{dollar}" "S{dollar}_1{dollar}" "S{dollar}_0{dollar} (EXACT)" "S{dollar}_0{dollar} (FSSH)" "S{dollar}_0{dollar} (LZSH)" "S{dollar}_0{dollar} (MASH)" --xlabel "Coordinate (a.u.)" "Time (a.u.)" --ylabel "Energy (a.u.)" "Ground State Population" --output SH_COMPARE.png')

# DYNAMICS ON URACIL VC MODEL ================================================================================================================================================================

MAITRA_MCTDH_D0_X = [0.00000000000, 32.40831373777, 87.25315237093, 149.57683263588, 191.95693521605, 234.33703779621, 289.18187642937, 334.05492622013, 353.99850390492, 383.91387043210, 413.82923695927, 456.20933953944, 486.12470606662, 513.54712538319, 530.99775585738, 560.91312238456, 585.84259449054, 613.26501380712, 675.58869407207, 732.92647991582, 785.27837133838, 880.01036534111, 952.30583444845, 999.67183144981, 1074.46024776776, 1134.29098082211, 1181.65697782348, 1243.98065808843, 1348.68444093355, 1438.43054051507, 1530.66958730720, 1607.95095083574, 1677.75347273249, 1752.54188905043, 1832.31619978957, 1897.13282726512, 2001.83661011023, 2061.66734316459, 2141.44165390372, 2193.79354532628, 2261.10312001243, 2325.91974748798, 2398.21521659532, 2477.98952733446]
MAITRA_MCTDH_D0_Y = [0.00911213548,  0.00998573466,  0.01997146932,   0.02710413694,   0.03708987161,   0.05848787446,   0.07417974322,   0.09415121255,   0.11840228245,   0.14835948644,   0.17403708987,   0.20542082738,   0.23395149786,   0.27246790299,   0.29671897289,   0.33380884450,   0.36376604850,   0.38373751783,   0.41084165477,   0.42938659058,   0.46932952924,   0.46504992867,   0.48074179743,   0.50499286733,    0.51925820256,    0.53780313837,    0.56062767475,    0.56776034236,    0.58059914407,    0.58915834522,    0.59486447931,    0.61198288159,    0.61198288159,    0.61198288159,    0.60485021398,    0.60485021398,    0.63052781740,    0.63195435092,    0.64907275320,    0.66333808844,    0.65905848787,    0.68616262482,    0.69900142653,    0.70613409415]

if URACIL_LVC:

    omega = {
        8  : np.array([734, 771, 1193, 1383, 1406, 1673, 1761, 1794]) / 8065.54429 / 27.211324570273
    }

    ap_input_8_10 = {"model_potential" : {"adiabatic" : True, "limits" : [-2.0, 2.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q10.mat", "points" : 4096, "potential" : "uracilDimless8D_1", "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 0]}}
    ap_input_8_12 = {"model_potential" : {"adiabatic" : True, "limits" : [-2.0, 2.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q12.mat", "points" : 4096, "potential" : "uracilDimless8D_1", "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 1]}}
    ap_input_8_18 = {"model_potential" : {"adiabatic" : True, "limits" : [-2.0, 2.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q18.mat", "points" : 4096, "potential" : "uracilDimless8D_1", "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 2]}}
    ap_input_8_20 = {"model_potential" : {"adiabatic" : True, "limits" : [-3.5, 3.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q20.mat", "points" : 4096, "potential" : "uracilDimless8D_1", "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 3]}}
    ap_input_8_21 = {"model_potential" : {"adiabatic" : True, "limits" : [-4.5, 4.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q21.mat", "points" : 4096, "potential" : "uracilDimless8D_1", "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 4]}}
    ap_input_8_24 = {"model_potential" : {"adiabatic" : True, "limits" : [-3.5, 3.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q24.mat", "points" : 4096, "potential" : "uracilDimless8D_1", "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 5]}}
    ap_input_8_25 = {"model_potential" : {"adiabatic" : True, "limits" : [-4.5, 4.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q25.mat", "points" : 4096, "potential" : "uracilDimless8D_1", "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 6]}}
    ap_input_8_26 = {"model_potential" : {"adiabatic" : True, "limits" : [-4.5, 4.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q26.mat", "points" : 4096, "potential" : "uracilDimless8D_1", "constant" : [{"index" : i, "value" : 0} for i in range(8) if i != 7]}}

    ap_input_12_3  = {"model_potential" : {"adiabatic" : True, "limits" : [-4.5, 4.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q3.mat",  "points" : 4096, "potential" : "uracilDimless12D_1", "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  0]}}
    ap_input_12_7  = {"model_potential" : {"adiabatic" : True, "limits" : [-3.0, 3.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q7.mat",  "points" : 4096, "potential" : "uracilDimless12D_1", "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  1]}}
    ap_input_12_10 = {"model_potential" : {"adiabatic" : True, "limits" : [-2.0, 2.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q10.mat", "points" : 4096, "potential" : "uracilDimless12D_1", "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  2]}}
    ap_input_12_11 = {"model_potential" : {"adiabatic" : True, "limits" : [-3.0, 3.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q11.mat", "points" : 4096, "potential" : "uracilDimless12D_1", "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  3]}}
    ap_input_12_12 = {"model_potential" : {"adiabatic" : True, "limits" : [-2.0, 2.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q12.mat", "points" : 4096, "potential" : "uracilDimless12D_1", "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  4]}}
    ap_input_12_18 = {"model_potential" : {"adiabatic" : True, "limits" : [-2.0, 2.0], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q18.mat", "points" : 4096, "potential" : "uracilDimless12D_1", "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  5]}}
    ap_input_12_19 = {"model_potential" : {"adiabatic" : True, "limits" : [-3.5, 3.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q19.mat", "points" : 4096, "potential" : "uracilDimless12D_1", "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  6]}}
    ap_input_12_20 = {"model_potential" : {"adiabatic" : True, "limits" : [-3.5, 3.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q20.mat", "points" : 4096, "potential" : "uracilDimless12D_1", "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  7]}}
    ap_input_12_21 = {"model_potential" : {"adiabatic" : True, "limits" : [-4.5, 4.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q21.mat", "points" : 4096, "potential" : "uracilDimless12D_1", "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  8]}}
    ap_input_12_24 = {"model_potential" : {"adiabatic" : True, "limits" : [-3.5, 3.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q24.mat", "points" : 4096, "potential" : "uracilDimless12D_1", "constant" : [{"index" : i, "value" : 0} for i in range(12) if i !=  9]}}
    ap_input_12_25 = {"model_potential" : {"adiabatic" : True, "limits" : [-4.5, 4.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q25.mat", "points" : 4096, "potential" : "uracilDimless12D_1", "constant" : [{"index" : i, "value" : 0} for i in range(12) if i != 10]}}
    ap_input_12_26 = {"model_potential" : {"adiabatic" : True, "limits" : [-4.5, 4.5], "output" : "URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q26.mat", "points" : 4096, "potential" : "uracilDimless12D_1", "constant" : [{"index" : i, "value" : 0} for i in range(12) if i != 11]}}

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
                "trajectories" : TRAJS
            }
        }

        fs_input = cp.deepcopy(cd_input); fs_input["classical_dynamics"]["fewest_switches"] = {}
        kt_input = cp.deepcopy(cd_input); kt_input["classical_dynamics"]["fewest_switches"] = {}
        lz_input = cp.deepcopy(cd_input); lz_input["classical_dynamics"]["landau_zener"]    = {}

        kt_input["classical_dynamics"]["time_derivative_coupling"] = "baeckan"

        fs_input["classical_dynamics"]["write"]["population_mean"] = f"URACIL_LVC_POPULATION_uracil{i}D_1_FSSH.mat"
        kt_input["classical_dynamics"]["write"]["population_mean"] = f"URACIL_LVC_POPULATION_uracil{i}D_1_KTSH.mat"
        lz_input["classical_dynamics"]["write"]["population_mean"] = f"URACIL_LVC_POPULATION_uracil{i}D_1_LZSH.mat"

        open(f"input_URACIL_LVC_fs_uracil8D_1.json", "w").write(js.dumps(fs_input, indent=4)); os.system(f"acorn input_URACIL_LVC_fs_uracil8D_1.json")
        open(f"input_URACIL_LVC_kt_uracil8D_1.json", "w").write(js.dumps(kt_input, indent=4)); os.system(f"acorn input_URACIL_LVC_kt_uracil8D_1.json")
        open(f"input_URACIL_LVC_lz_uracil8D_1.json", "w").write(js.dumps(lz_input, indent=4)); os.system(f"acorn input_URACIL_LVC_lz_uracil8D_1.json")

        np.savetxt("URACIL_LVC_MAITRA_MCTDH_D0.mat", np.array([MAITRA_MCTDH_D0_X, MAITRA_MCTDH_D0_Y]).T, comments="", fmt="%20.14f", header=f"{len(MAITRA_MCTDH_D0_X)} 2")

        os.system(f'python python/lines.py URACIL_LVC_MAITRA_MCTDH_D0.mat URACIL_LVC_POPULATION_uracil{i}D_1_FSSH.mat:0 URACIL_LVC_POPULATION_uracil{i}D_1_KTSH.mat:0 URACIL_LVC_POPULATION_uracil{i}D_1_LZSH.mat:0 --legends every "D{dollar}_0{dollar} (MCTDH)" "D{dollar}_0{dollar} (FSSH)" "D{dollar}_0{dollar} ({dollar}{esc}kappa{dollar}TSH)" "D{dollar}_0{dollar} (LZSH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_uracil{i}D_1.png')

    os.system(f'python python/lines.py URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q18.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q20.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q21.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q24.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q25.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q26.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q10.mat:0,5,10 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless8D_1_Q12.mat:0,5,10 --figsize 6 8 --subplots 241 241 241 241 242 242 242 242 243 243 243 243 244 244 244 244 245 245 245 245 246 246 246 246 247 247 247 248 248 248 --xlabel "Q{dollar}_{{18}}{dollar}" "Q{dollar}_{{20}}{dollar}" "Q{dollar}_{{21}}{dollar}" "Q{dollar}_{{24}}{dollar}" "Q{dollar}_{{25}}{dollar}" "Q{dollar}_{{26}}{dollar}" "Q{dollar}_{{10}}{dollar}" "Q{dollar}_{{12}}{dollar}" --output POTENTIAL_ADIABATIC_uracilDimless8D_1.png')
    os.system(f'python python/lines.py URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q3.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q7.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q11.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q18.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q19.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q20.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q21.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q24.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q25.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q26.mat:0,5,10,15 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q10.mat:0,5,10 URACIL_LVC_POTENTIAL_ADIABATIC_uracilDimless12D_1_Q12.mat:0,5,10 --figsize 9 8 --subplots 341 341 341 341 342 342 342 342 343 343 343 343 344 344 344 344 345 345 345 345 346 346 346 346 347 347 347 347 348 348 348 348 349 349 349 349 3-4-10 3-4-10 3-4-10 3-4-10 3-4-11 3-4-11 3-4-11 3-4-12 3-4-12 3-4-12 --xlabel "Q{dollar}_{{3}}{dollar}" "Q{dollar}_{{7}}{dollar}" "Q{dollar}_{{11}}{dollar}" "Q{dollar}_{{18}}{dollar}" "Q{dollar}_{{19}}{dollar}" "Q{dollar}_{{20}}{dollar}" "Q{dollar}_{{21}}{dollar}" "Q{dollar}_{{24}}{dollar}" "Q{dollar}_{{25}}{dollar}" "Q{dollar}_{{26}}{dollar}" "Q{dollar}_{{10}}{dollar}" "Q{dollar}_{{12}}{dollar}" --output POTENTIAL_ADIABATIC_uracilDimless12D_1.png')
# CLEANUP ====================================================================================================================================================================================

if os.path.exists("research"): sh.rmtree("research")

os.mkdir("research"); os.mkdir("research/source"); os.mkdir("research/output")

for file in gl.glob("*.json"): sh.move(file, "research/source")
for file in gl.glob("*.mat" ): sh.move(file, "research/source")
for file in gl.glob("*.png" ): sh.move(file, "research/output")
