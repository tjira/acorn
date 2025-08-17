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

        .iterations = 1000,
        .time_step = 0.01,
        .trajectories = 1000,
    };

    const output_harmonic1D_1 = try classical_dynamics.run(f64, opt_harmonic1D_1, false, allocator); defer output_harmonic1D_1.deinit();

    try expect(output_harmonic1D_1.r.at(0),                         -0.86782025746194);
    try expect(output_harmonic1D_1.p.at(0),                          0.51552062843856);
    try expect(output_harmonic1D_1.Ekin + output_harmonic1D_1.Epot,  1.10923079930487);
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

        .iterations = 3500,
        .time_step = 1,
        .trajectories = 100,
    };

    var opt_tully1D_2 = opt_tully1D_1; opt_tully1D_2.hamiltonian.name = "tully1D_2";
    var opt_tully1D_3 = opt_tully1D_1; opt_tully1D_3.hamiltonian.name = "tully1D_3";

    const output_tully1D_1 = try classical_dynamics.run(f64, opt_tully1D_1, false, allocator); defer output_tully1D_1.deinit();
    const output_tully1D_2 = try classical_dynamics.run(f64, opt_tully1D_2, false, allocator); defer output_tully1D_2.deinit();
    const output_tully1D_3 = try classical_dynamics.run(f64, opt_tully1D_3, false, allocator); defer output_tully1D_3.deinit();

    try expect(output_tully1D_1.r.at(0)                     , 20.94084643908319);
    try expect(output_tully1D_1.p.at(0)                     , 20.99332774403862);
    try expect(output_tully1D_1.Ekin + output_tully1D_1.Epot,  0.10955586987758);
    try expect(output_tully1D_1.pop.at(0)                   ,  0.55000000000000);
    try expect(output_tully1D_1.pop.at(1)                   ,  0.45000000000000);

    try expect(output_tully1D_2.r.at(0)                     , 21.39647898664311);
    try expect(output_tully1D_2.p.at(0)                     , 20.64948572850306);
    try expect(output_tully1D_2.Ekin + output_tully1D_2.Epot,  0.14962938883238);
    try expect(output_tully1D_2.pop.at(0)                   ,  0.16000000000000);
    try expect(output_tully1D_2.pop.at(1)                   ,  0.84000000000000);

    try expect(output_tully1D_3.r.at(0)                     ,  1.39348081243005);
    try expect(output_tully1D_3.p.at(0)                     , -0.79439349967970);
    try expect(output_tully1D_3.Ekin + output_tully1D_3.Epot,  0.10015197887868);
    try expect(output_tully1D_3.pop.at(0)                   ,  0.58000000000000);
    try expect(output_tully1D_3.pop.at(1)                   ,  0.42000000000000);
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

        .iterations = 3500,
        .time_step = 1,
        .trajectories = 100,
    };

    var opt_tully1D_2 = opt_tully1D_1; opt_tully1D_2.hamiltonian.name = "tully1D_2";
    var opt_tully1D_3 = opt_tully1D_1; opt_tully1D_3.hamiltonian.name = "tully1D_3";

    const output_tully1D_1 = try classical_dynamics.run(f64, opt_tully1D_1, false, allocator); defer output_tully1D_1.deinit();
    const output_tully1D_2 = try classical_dynamics.run(f64, opt_tully1D_2, false, allocator); defer output_tully1D_2.deinit();
    const output_tully1D_3 = try classical_dynamics.run(f64, opt_tully1D_3, false, allocator); defer output_tully1D_3.deinit();

    try expect(output_tully1D_1.r.at(0)                     ,  21.03388137082108);
    try expect(output_tully1D_1.p.at(0)                     ,  21.11075780135698);
    try expect(output_tully1D_1.Ekin + output_tully1D_1.Epot,   0.10954971275913);
    try expect(output_tully1D_1.pop.at(0)                   ,   0.61000000000000);
    try expect(output_tully1D_1.pop.at(1)                   ,   0.39000000000000);

    try expect(output_tully1D_2.r.at(0)                     ,  21.35474574219164);
    try expect(output_tully1D_2.p.at(0)                     ,  20.66377958780642);
    try expect(output_tully1D_2.Ekin + output_tully1D_2.Epot,   0.14955233291760);
    try expect(output_tully1D_2.pop.at(0)                   ,   0.16000000000000);
    try expect(output_tully1D_2.pop.at(1)                   ,   0.84000000000000);

    try expect(output_tully1D_3.r.at(0)                     , -16.75966919104387);
    try expect(output_tully1D_3.p.at(0)                     , -19.93694370568629);
    try expect(output_tully1D_3.Ekin + output_tully1D_3.Epot,   0.10015338907023);
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

        .iterations = 3500,
        .time_step = 1,
        .trajectories = 100,
    };

    var opt_tully1D_2 = opt_tully1D_1; opt_tully1D_2.hamiltonian.name = "tully1D_2";
    var opt_tully1D_3 = opt_tully1D_1; opt_tully1D_3.hamiltonian.name = "tully1D_3";

    const output_tully1D_1 = try classical_dynamics.run(f64, opt_tully1D_1, false, allocator); defer output_tully1D_1.deinit();
    const output_tully1D_2 = try classical_dynamics.run(f64, opt_tully1D_2, false, allocator); defer output_tully1D_2.deinit();
    const output_tully1D_3 = try classical_dynamics.run(f64, opt_tully1D_3, false, allocator); defer output_tully1D_3.deinit();

    try expect(output_tully1D_1.r.at(0)                     , 20.84372569926736);
    try expect(output_tully1D_1.p.at(0)                     , 20.89958627219347);
    try expect(output_tully1D_1.Ekin + output_tully1D_1.Epot,  0.10955004442831);
    try expect(output_tully1D_1.pop.at(0)                   ,  0.50000000000000);
    try expect(output_tully1D_1.pop.at(1)                   ,  0.50000000000000);

    try expect(output_tully1D_2.r.at(0)                     , 21.95063987602763);
    try expect(output_tully1D_2.p.at(0)                     , 21.22225043142191);
    try expect(output_tully1D_2.Ekin + output_tully1D_2.Epot,  0.14961662501432);
    try expect(output_tully1D_2.pop.at(0)                   ,  0.28000000000000);
    try expect(output_tully1D_2.pop.at(1)                   ,  0.72000000000000);

    try expect(output_tully1D_3.r.at(0)                     ,  9.45525584238022);
    try expect(output_tully1D_3.p.at(0)                     ,  7.90041858544399);
    try expect(output_tully1D_3.Ekin + output_tully1D_3.Epot,  0.10015169692535);
    try expect(output_tully1D_3.pop.at(0)                   ,  0.71000000000000);
    try expect(output_tully1D_3.pop.at(1)                   ,  0.29000000000000);
}
