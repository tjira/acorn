const std = @import("std");

const input              = @import("acorn").input             ;
const classical_dynamics = @import("acorn").classical_dynamics;

const expect = @import("main.zig").expect;
const log    = @import("main.zig").log   ;

const allocator = std.testing.allocator;

test "classical_dynamics_adiabatic" {
    const opt_harmonic1D_1 = input.ClassicalDynamicsOptions(f64){
        .hamiltonian = .{
            .name = "harmonic1D_1"
        },
        .initial_conditions = .{
            .position_mean = &[_]f64{1.0},
            .position_std  = &[_]f64{0.5},
            .momentum_mean = &[_]f64{0.0},
            .momentum_std  = &[_]f64{1.0},
            .state = &[_]f64{1}, .mass = &[_]f64{1}
        },

        .adiabatic = true,
        .iterations = 1000,
        .time_step = 0.01,
        .trajectories = 1000,
    };

    const output_harmonic1D_1 = try classical_dynamics.run(f64, opt_harmonic1D_1, false, allocator); defer output_harmonic1D_1.deinit();

    try expect(output_harmonic1D_1.r.at(0),                         -0.86782025746194);
    try expect(output_harmonic1D_1.p.at(0),                          0.51552062843856);
    try expect(output_harmonic1D_1.Ekin + output_harmonic1D_1.Epot,  1.10836247904741);
}

test "classical_dynamics_nonadiabatic_fssh" {
    const opt_tully1D_1 = input.ClassicalDynamicsOptions(f64){
        .hamiltonian = .{
            .name = "tully1D_1"
        },
        .initial_conditions = .{
            .position_mean = &[_]f64{-15.0},
            .position_std  = &[_]f64{  0.5},
            .momentum_mean = &[_]f64{ 20.0},
            .momentum_std  = &[_]f64{  1.0},
            .state = &[_]f64{0, 1}, .mass = &[_]f64{2000}
        },

        .fewest_switches = .{},

        .adiabatic = true,
        .iterations = 3500,
        .time_step = 1,
        .trajectories = 100,
    };

    var opt_tully1D_2 = opt_tully1D_1; opt_tully1D_2.hamiltonian.name = "tully1D_2";
    var opt_tully1D_3 = opt_tully1D_1; opt_tully1D_3.hamiltonian.name = "tully1D_3";

    const output_tully1D_1 = try classical_dynamics.run(f64, opt_tully1D_1, false, allocator); defer output_tully1D_1.deinit();
    const output_tully1D_2 = try classical_dynamics.run(f64, opt_tully1D_2, false, allocator); defer output_tully1D_2.deinit();
    const output_tully1D_3 = try classical_dynamics.run(f64, opt_tully1D_3, false, allocator); defer output_tully1D_3.deinit();

    try expect(output_tully1D_1.r.at(0)                     , 20.94080629806638);
    try expect(output_tully1D_1.p.at(0)                     , 20.99328175367489);
    try expect(output_tully1D_1.Ekin + output_tully1D_1.Epot,  0.10955539972974);
    try expect(output_tully1D_1.pop.at(0)                   ,  0.55000000000000);
    try expect(output_tully1D_1.pop.at(1)                   ,  0.45000000000000);

    try expect(output_tully1D_2.r.at(0)                     , 21.39557352463867);
    try expect(output_tully1D_2.p.at(0)                     , 20.64848889542154);
    try expect(output_tully1D_2.Ekin + output_tully1D_2.Epot,  0.14961909128636);
    try expect(output_tully1D_2.pop.at(0)                   ,  0.16000000000000);
    try expect(output_tully1D_2.pop.at(1)                   ,  0.84000000000000);

    try expect(output_tully1D_3.r.at(0)                     ,  5.47747352257774);
    try expect(output_tully1D_3.p.at(0)                     ,  3.54994291467015);
    try expect(output_tully1D_3.Ekin + output_tully1D_3.Epot,  0.10015363583858);
    try expect(output_tully1D_3.pop.at(0)                   ,  0.63000000000000);
    try expect(output_tully1D_3.pop.at(1)                   ,  0.37000000000000);
}

test "classical_dynamics_nonadiabatic_lzsh" {
    const opt_tully1D_1 = input.ClassicalDynamicsOptions(f64){
        .hamiltonian = .{
            .name = "tully1D_1"
        },
        .initial_conditions = .{
            .position_mean = &[_]f64{-15.0},
            .position_std  = &[_]f64{  0.5},
            .momentum_mean = &[_]f64{ 20.0},
            .momentum_std  = &[_]f64{  1.0},
            .state = &[_]f64{0, 1}, .mass = &[_]f64{2000}
        },

        .landau_zener = .{},

        .adiabatic = true,
        .iterations = 3500,
        .time_step = 1,
        .trajectories = 100,
    };

    var opt_tully1D_2 = opt_tully1D_1; opt_tully1D_2.hamiltonian.name = "tully1D_2";
    var opt_tully1D_3 = opt_tully1D_1; opt_tully1D_3.hamiltonian.name = "tully1D_3";

    const output_tully1D_1 = try classical_dynamics.run(f64, opt_tully1D_1, false, allocator); defer output_tully1D_1.deinit();
    const output_tully1D_2 = try classical_dynamics.run(f64, opt_tully1D_2, false, allocator); defer output_tully1D_2.deinit();
    const output_tully1D_3 = try classical_dynamics.run(f64, opt_tully1D_3, false, allocator); defer output_tully1D_3.deinit();

    try expect(output_tully1D_1.r.at(0)                     ,  21.03394367044433);
    try expect(output_tully1D_1.p.at(0)                     ,  21.11081086713970);
    try expect(output_tully1D_1.Ekin + output_tully1D_1.Epot,   0.10955029365272);
    try expect(output_tully1D_1.pop.at(0)                   ,   0.61000000000000);
    try expect(output_tully1D_1.pop.at(1)                   ,   0.39000000000000);

    try expect(output_tully1D_2.r.at(0)                     ,  21.35475558355247);
    try expect(output_tully1D_2.p.at(0)                     ,  20.66378942157508);
    try expect(output_tully1D_2.Ekin + output_tully1D_2.Epot,   0.14955245323037);
    try expect(output_tully1D_2.pop.at(0)                   ,   0.16000000000000);
    try expect(output_tully1D_2.pop.at(1)                   ,   0.84000000000000);

    try expect(output_tully1D_3.r.at(0)                     , -16.75966919104387);
    try expect(output_tully1D_3.p.at(0)                     , -19.93694370568629);
    try expect(output_tully1D_3.Ekin + output_tully1D_3.Epot,   0.10015338907028);
    try expect(output_tully1D_3.pop.at(0)                   ,   0.00000000000000);
    try expect(output_tully1D_3.pop.at(1)                   ,   1.00000000000000);
}

test "classical_dynamics_nonadiabatic_mash" {
    const opt_tully1D_1 = input.ClassicalDynamicsOptions(f64){
        .hamiltonian = .{
            .name = "tully1D_1"
        },
        .initial_conditions = .{
            .position_mean = &[_]f64{-15.0},
            .position_std  = &[_]f64{  0.5},
            .momentum_mean = &[_]f64{ 20.0},
            .momentum_std  = &[_]f64{  1.0},
            .state = &[_]f64{0, 1}, .mass = &[_]f64{2000}
        },

        .spin_mapping = .{},

        .adiabatic = true,
        .iterations = 3500,
        .time_step = 1,
        .trajectories = 100,
    };

    var opt_tully1D_2 = opt_tully1D_1; opt_tully1D_2.hamiltonian.name = "tully1D_2";
    var opt_tully1D_3 = opt_tully1D_1; opt_tully1D_3.hamiltonian.name = "tully1D_3";

    const output_tully1D_1 = try classical_dynamics.run(f64, opt_tully1D_1, false, allocator); defer output_tully1D_1.deinit();
    const output_tully1D_2 = try classical_dynamics.run(f64, opt_tully1D_2, false, allocator); defer output_tully1D_2.deinit();
    const output_tully1D_3 = try classical_dynamics.run(f64, opt_tully1D_3, false, allocator); defer output_tully1D_3.deinit();

    try expect(output_tully1D_1.r.at(0)                     , 21.00843552107591);
    try expect(output_tully1D_1.p.at(0)                     , 21.06108731339015);
    try expect(output_tully1D_1.Ekin + output_tully1D_1.Epot,  0.10955962058843);
    try expect(output_tully1D_1.pop.at(0)                   ,  0.59000000000000);
    try expect(output_tully1D_1.pop.at(1)                   ,  0.41000000000000);

    try expect(output_tully1D_2.r.at(0)                     , 21.11749313656524);
    try expect(output_tully1D_2.p.at(0)                     , 20.38803787851100);
    try expect(output_tully1D_2.Ekin + output_tully1D_2.Epot,  0.14959614430443);
    try expect(output_tully1D_2.pop.at(0)                   ,  0.10000000000000);
    try expect(output_tully1D_2.pop.at(1)                   ,  0.90000000000000);

    try expect(output_tully1D_3.r.at(0)                     , -1.14435520049919);
    try expect(output_tully1D_3.p.at(0)                     , -3.60940500235558);
    try expect(output_tully1D_3.Ekin + output_tully1D_3.Epot,  0.10015587901201);
    try expect(output_tully1D_3.pop.at(0)                   ,  0.68000000000000);
    try expect(output_tully1D_3.pop.at(1)                   ,  0.32000000000000);
}
