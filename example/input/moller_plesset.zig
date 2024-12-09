{
    "moller_plesset" : {
        "hartree_fock" : {
            "molecule" : "example/molecule/water.xyz",
            "integral" : {
                "overlap" : "S_AO.mat",
                "kinetic" : "T_AO.mat",
                "nuclear" : "V_AO.mat",
                "coulomb" : "J_AO.mat"
            },
            "threshold" : 1e-8,
            "maxiter" : 100
        },
        "order" : 2
    }
}
