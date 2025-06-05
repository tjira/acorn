const std = @import("std");

const input            = @import("acorn").input;
const quantum_dynamics = @import("acorn").quantum_dynamics;

const expect = @import("main.zig").expect;
const log    = @import("main.zig").log   ;

const allocator = std.testing.allocator;

test "quantum_dynamics_imaginary_adiabatic-1d" {
    const opt_harmonic1D_1 = input.QuantumDynamicsOptions(f64){
        .grid = .{
            .limits = &[_]f64{-8, 8}, .points = 64
        },
        .hamiltonian = .{
            .name = "harmonic1D_1"
        },
        .initial_conditions = .{
            .position = &[_]f64{1}, .momentum = &[_]f64{0}, .gamma = 2, .state = 0, .mass = 1
        },

        .adiabatic = false,
        .iterations = 1000,
        .mode = &[_]u32{10, 0},
        .time_step = 0.01,
    };

    const output_harmonic1D_1 = try quantum_dynamics.run(f64, opt_harmonic1D_1, false, allocator); defer output_harmonic1D_1.deinit();

    try expect(output_harmonic1D_1.r[0].at(0),  0.00006053550963);
    try expect(output_harmonic1D_1.r[1].at(0), -0.00002269862596);
    try expect(output_harmonic1D_1.r[2].at(0), -0.00004994211385);
    try expect(output_harmonic1D_1.r[3].at(0),  0.00041323220537);
    try expect(output_harmonic1D_1.r[4].at(0), -0.00035429312937);
    try expect(output_harmonic1D_1.r[5].at(0), -0.00013301655239);
    try expect(output_harmonic1D_1.r[6].at(0),  0.00024239920495);
    try expect(output_harmonic1D_1.r[7].at(0), -0.00015726015675);
    try expect(output_harmonic1D_1.r[8].at(0),  0.01059601015868);
    try expect(output_harmonic1D_1.r[9].at(0), -0.01053559862599);

    try expect(output_harmonic1D_1.p[0].at(0), -0.00000000000003);
    try expect(output_harmonic1D_1.p[1].at(0), -0.00000000000001);
    try expect(output_harmonic1D_1.p[2].at(0),  0.00000000000001);
    try expect(output_harmonic1D_1.p[3].at(0),  0.00000000000000);
    try expect(output_harmonic1D_1.p[4].at(0), -0.00000000000001);
    try expect(output_harmonic1D_1.p[5].at(0), -0.00000000000000);
    try expect(output_harmonic1D_1.p[6].at(0),  0.00000000000001);
    try expect(output_harmonic1D_1.p[7].at(0),  0.00000000000000);
    try expect(output_harmonic1D_1.p[8].at(0),  0.00000000000001);
    try expect(output_harmonic1D_1.p[9].at(0), -0.00000000000004);

    try expect(output_harmonic1D_1.Ekin[0] + output_harmonic1D_1.Epot[0], 0.50000000187134);
    try expect(output_harmonic1D_1.Ekin[1] + output_harmonic1D_1.Epot[1], 1.49999999864282);
    try expect(output_harmonic1D_1.Ekin[2] + output_harmonic1D_1.Epot[2], 2.49999999986185);
    try expect(output_harmonic1D_1.Ekin[3] + output_harmonic1D_1.Epot[3], 3.50000002036197);
    try expect(output_harmonic1D_1.Ekin[4] + output_harmonic1D_1.Epot[4], 4.49999998045793);
    try expect(output_harmonic1D_1.Ekin[5] + output_harmonic1D_1.Epot[5], 5.50000000082936);
    try expect(output_harmonic1D_1.Ekin[6] + output_harmonic1D_1.Epot[6], 6.50000000163189);
    try expect(output_harmonic1D_1.Ekin[7] + output_harmonic1D_1.Epot[7], 7.49999999884292);
    try expect(output_harmonic1D_1.Ekin[8] + output_harmonic1D_1.Epot[8], 8.50000623707168);
    try expect(output_harmonic1D_1.Ekin[9] + output_harmonic1D_1.Epot[9], 9.49999376451070);
}

test "quantum_dynamics_real_adiabatic-1d" {
    const opt_harmonic1D_1 = input.QuantumDynamicsOptions(f64){
        .grid = .{
            .limits = &[_]f64{-8, 8}, .points = 64
        },
        .hamiltonian = .{
            .name = "harmonic1D_1"
        },
        .initial_conditions = .{
            .position = &[_]f64{1}, .momentum = &[_]f64{0}, .gamma = 2, .state = 0, .mass = 1
        },

        .adiabatic = false,
        .iterations = 1000,
        .mode = &[_]u32{0, 1},
        .time_step = 0.01,
    };

    const output_harmonic1D_1 = try quantum_dynamics.run(f64, opt_harmonic1D_1, false, allocator); defer output_harmonic1D_1.deinit();

    try expect(output_harmonic1D_1.r[0].at(0)                           , -0.83904886054760);
    try expect(output_harmonic1D_1.p[0].at(0)                           ,  0.54404927138121);
    try expect(output_harmonic1D_1.Ekin[0] + output_harmonic1D_1.Epot[0],  1.12499907510295);
}

test "quantum_dynamics_real_nonadiabatic-1d/2s" {
    const opt_tully1D_1 = input.QuantumDynamicsOptions(f64){
        .grid = .{
            .limits = &[_]f64{-16, 32}, .points = 2048
        },
        .hamiltonian = .{
            .name = "tully1D_1"
        },
        .initial_conditions = .{
            .position = &[_]f64{-15}, .momentum = &[_]f64{20}, .gamma = 2, .state = 1, .mass = 2000
        },

        .adiabatic = true,
        .iterations = 350,
        .mode = &[_]u32{0, 1},
        .time_step = 10,
    };

    var opt_tully1D_2 = opt_tully1D_1; opt_tully1D_2.hamiltonian.name = "tully1D_2";
    var opt_tully1D_3 = opt_tully1D_1; opt_tully1D_3.hamiltonian.name = "tully1D_3";

    const output_tully1D_1 = try quantum_dynamics.run(f64, opt_tully1D_1, false, allocator); defer output_tully1D_1.deinit();
    const output_tully1D_2 = try quantum_dynamics.run(f64, opt_tully1D_2, false, allocator); defer output_tully1D_2.deinit();
    const output_tully1D_3 = try quantum_dynamics.run(f64, opt_tully1D_3, false, allocator); defer output_tully1D_3.deinit();

    try expect(output_tully1D_1.r[0].at(0)                        , 21.01938968491557);
    try expect(output_tully1D_1.p[0].at(0)                        , 21.02856421418769);
    try expect(output_tully1D_1.Ekin[0] + output_tully1D_1.Epot[0],  0.11177076998784);
    try expect(output_tully1D_1.P[0].at(0, 0).re                  ,  0.54705439670020);
    try expect(output_tully1D_1.P[0].at(0, 1).re                  ,  0.09068606015860);
    try expect(output_tully1D_1.P[0].at(1, 0).re                  ,  0.09068606015860);
    try expect(output_tully1D_1.P[0].at(1, 1).re                  ,  0.45294560329986);

    try expect(output_tully1D_2.r[0].at(0)                        , 21.58373845982900);
    try expect(output_tully1D_2.p[0].at(0)                        , 20.83738069203519);
    try expect(output_tully1D_2.Ekin[0] + output_tully1D_2.Epot[0],  0.15176981510263);
    try expect(output_tully1D_2.P[0].at(0, 0).re                  ,  0.19319230897050);
    try expect(output_tully1D_2.P[0].at(0, 1).re                  , -0.00121002615618);
    try expect(output_tully1D_2.P[0].at(1, 0).re                  , -0.00121002615618);
    try expect(output_tully1D_2.P[0].at(1, 1).re                  ,  0.80680769102958);

    try expect(output_tully1D_3.r[0].at(0)                        , -3.63335060067555);
    try expect(output_tully1D_3.p[0].at(0)                        ,  9.79711483366324);
    try expect(output_tully1D_3.Ekin[0] + output_tully1D_3.Epot[0],  0.39648649056319);
    try expect(output_tully1D_3.P[0].at(0, 0).re                  ,  0.65976772649603);
    try expect(output_tully1D_3.P[0].at(0, 1).re                  , -0.18336140276161);
    try expect(output_tully1D_3.P[0].at(1, 0).re                  , -0.18336140276161);
    try expect(output_tully1D_3.P[0].at(1, 1).re                  ,  0.34023227350403);
}
