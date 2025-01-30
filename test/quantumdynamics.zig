const std = @import("std");

const quantum_dynamics = @import("acorn").quantum_dynamics;

const allocator = std.testing.allocator;

test "quantum_dynamics_imaginary_adiabatic-1d" {
    const opt_harmonic1D_1 = quantum_dynamics.QuantumDynamicsOptions(f64){
        .grid = .{
            .limits = &[_]f64{-8, 8}, .points = 64
        },
        .initial_conditions = .{
            .position = &[_]f64{1}, .momentum = &[_]f64{0}, .gamma = 2, .state = 0, .mass = 1
        },

        .adiabatic = false,
        .iterations = 1000,
        .mode = &[_]u32{10, 0},
        .potential = "harmonic1D_1",
        .time_step = 0.01,
    };

    const output_harmonic1D_1 = try quantum_dynamics.run(f64, opt_harmonic1D_1, false, allocator); defer output_harmonic1D_1.deinit();

    try std.testing.expect(@abs(output_harmonic1D_1.r[0].at(0) - 0.00006053550963) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.r[1].at(0) + 0.00002269862596) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.r[2].at(0) + 0.00004994211385) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.r[3].at(0) - 0.00041323220537) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.r[4].at(0) + 0.00035429312937) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.r[5].at(0) + 0.00013301655239) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.r[6].at(0) - 0.00024239920495) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.r[7].at(0) + 0.00015726015675) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.r[8].at(0) - 0.01059601015868) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.r[9].at(0) + 0.01053559862599) < 1e-14);

    try std.testing.expect(@abs(output_harmonic1D_1.p[0].at(0) + 0.00000000000003) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.p[1].at(0) + 0.00000000000001) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.p[2].at(0) - 0.00000000000001) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.p[3].at(0) - 0.00000000000000) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.p[4].at(0) + 0.00000000000001) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.p[5].at(0) + 0.00000000000001) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.p[6].at(0) - 0.00000000000001) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.p[7].at(0) - 0.00000000000000) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.p[8].at(0) + 0.00000000000005) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.p[9].at(0) - 0.00000000000002) < 1e-14);

    try std.testing.expect(@abs(output_harmonic1D_1.Ekin[0] + output_harmonic1D_1.Epot[0] - 0.50000000187134) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.Ekin[1] + output_harmonic1D_1.Epot[1] - 1.49999999864282) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.Ekin[2] + output_harmonic1D_1.Epot[2] - 2.49999999986185) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.Ekin[3] + output_harmonic1D_1.Epot[3] - 3.50000002036197) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.Ekin[4] + output_harmonic1D_1.Epot[4] - 4.49999998045793) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.Ekin[5] + output_harmonic1D_1.Epot[5] - 5.50000000082936) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.Ekin[6] + output_harmonic1D_1.Epot[6] - 6.50000000163189) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.Ekin[7] + output_harmonic1D_1.Epot[7] - 7.49999999884292) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.Ekin[8] + output_harmonic1D_1.Epot[8] - 8.50000623707168) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.Ekin[9] + output_harmonic1D_1.Epot[9] - 9.49999376451070) < 1e-14);
}

test "quantum_dynamics_real_adiabatic-1d" {
    const opt_harmonic1D_1 = quantum_dynamics.QuantumDynamicsOptions(f64){
        .grid = .{
            .limits = &[_]f64{-8, 8}, .points = 64
        },
        .initial_conditions = .{
            .position = &[_]f64{1}, .momentum = &[_]f64{0}, .gamma = 2, .state = 0, .mass = 1
        },

        .adiabatic = false,
        .iterations = 1000,
        .mode = &[_]u32{0, 1},
        .potential = "harmonic1D_1",
        .time_step = 0.01,
    };

    const output_harmonic1D_1 = try quantum_dynamics.run(f64, opt_harmonic1D_1, false, allocator); defer output_harmonic1D_1.deinit();

    try std.testing.expect(@abs(output_harmonic1D_1.r[0].at(0)                            + 0.83904886054760) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.p[0].at(0)                            - 0.54404927138121) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.Ekin[0] + output_harmonic1D_1.Epot[0] - 1.12499907510295) < 1e-14);
}

test "quantum_dynamics_real_nonadiabatic-1d/1s" {
    const opt_tully1D_1 = quantum_dynamics.QuantumDynamicsOptions(f64){
        .grid = .{
            .limits = &[_]f64{-16, 32}, .points = 2048
        },
        .initial_conditions = .{
            .position = &[_]f64{-10}, .momentum = &[_]f64{15}, .gamma = 2, .state = 1, .mass = 2000
        },

        .adiabatic = true,
        .iterations = 350,
        .mode = &[_]u32{0, 1},
        .potential = "tully1D_1",
        .time_step = 10,
    };

    var opt_tully1D_2 = opt_tully1D_1; opt_tully1D_2.potential = "tully1D_2"; opt_tully1D_2.initial_conditions.state = 1;
    var opt_tully1D_3 = opt_tully1D_1; opt_tully1D_3.potential = "tully1D_3"; opt_tully1D_3.initial_conditions.state = 0;

    const output_tully1D_1 = try quantum_dynamics.run(f64, opt_tully1D_1, false, allocator); defer output_tully1D_1.deinit();
    const output_tully1D_2 = try quantum_dynamics.run(f64, opt_tully1D_2, false, allocator); defer output_tully1D_2.deinit();
    const output_tully1D_3 = try quantum_dynamics.run(f64, opt_tully1D_3, false, allocator); defer output_tully1D_3.deinit();

    try std.testing.expect(@abs(output_tully1D_1.r[0].at(0)                         - 17.41407850665638) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_1.p[0].at(0)                         - 16.01050304057247) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_1.Ekin[0] + output_tully1D_1.Epot[0] -  0.06650000244591) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_1.P[0].at(0, 0)                      -  0.41050571809112) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_1.P[0].at(0, 1)                      +  0.03617449781705) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_1.P[0].at(1, 0)                      +  0.03617449781705) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_1.P[0].at(1, 1)                      -  0.58949428190631) < 1e-14);

    try std.testing.expect(@abs(output_tully1D_2.r[0].at(0)                         - 17.66897173759843) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_2.p[0].at(0)                         - 15.14618136068045) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_2.Ekin[0] + output_tully1D_2.Epot[0] -  0.10649999998912) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_2.P[0].at(0, 0)                      -  0.02483439803093) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_2.P[0].at(0, 1)                      -  0.00003143380659) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_2.P[0].at(1, 0)                      -  0.00003143380659) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_2.P[0].at(1, 1)                      -  0.97516560196651) < 1e-14);

    try std.testing.expect(@abs(output_tully1D_3.r[0].at(0)                         + 8.88235020923131) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_3.p[0].at(0)                         + 1.83824558248584) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_3.Ekin[0] + output_tully1D_3.Epot[0] - 0.23000238131947) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_3.P[0].at(0, 0)                      - 0.49287571352242) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_3.P[0].at(0, 1)                      + 0.13226570241379) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_3.P[0].at(1, 0)                      + 0.13226570241379) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_3.P[0].at(1, 1)                      - 0.50712428647487) < 1e-14);
}
