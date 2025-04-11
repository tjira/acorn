const std = @import("std");

const input              = @import("acorn").input             ;
const classical_dynamics = @import("acorn").classical_dynamics;

const expect = @import("main.zig").expect;
const log    = @import("main.zig").log   ;

const allocator = std.testing.allocator;

test "classical_dynamics_adiabatic" {
    const opt_harmonic1D_1 = input.ClassicalDynamicsOptions(f64){
        .initial_conditions = .{
            .position_mean = &[_]f64{1.0},
            .position_std  = &[_]f64{0.5},
            .momentum_mean = &[_]f64{0.0},
            .momentum_std  = &[_]f64{1.0},
            .state = &[_]f64{1}, .mass = &[_]f64{1}
        },

        .adiabatic = true,
        .iterations = 1000,
        .potential = "harmonic1D_1",
        .time_step = 0.01,
        .trajectories = 1000,
    };

    const output_harmonic1D_1 = try classical_dynamics.run(f64, opt_harmonic1D_1, false, allocator); defer output_harmonic1D_1.deinit();

    try expect(output_harmonic1D_1.r.at(0),                         -0.86233273468923);
    try expect(output_harmonic1D_1.p.at(0),                          0.51533412567938);
    try expect(output_harmonic1D_1.Ekin + output_harmonic1D_1.Epot,  1.10694105966608);
}

test "classical_dynamics_nonadiabatic_fssh" {
    const opt_tully1D_1 = input.ClassicalDynamicsOptions(f64){
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
        .potential = "tully1D_1",
        .time_step = 1,
        .trajectories = 100,
    };

    var opt_tully1D_2 = opt_tully1D_1; opt_tully1D_2.potential = "tully1D_2";
    var opt_tully1D_3 = opt_tully1D_1; opt_tully1D_3.potential = "tully1D_3";

    const output_tully1D_1 = try classical_dynamics.run(f64, opt_tully1D_1, false, allocator); defer output_tully1D_1.deinit();
    const output_tully1D_2 = try classical_dynamics.run(f64, opt_tully1D_2, false, allocator); defer output_tully1D_2.deinit();
    const output_tully1D_3 = try classical_dynamics.run(f64, opt_tully1D_3, false, allocator); defer output_tully1D_3.deinit();

    try expect(output_tully1D_1.r.at(0)                     , 20.96747749693828);
    try expect(output_tully1D_1.p.at(0)                     , 21.01999888560924);
    try expect(output_tully1D_1.Ekin + output_tully1D_1.Epot,  0.10954931488450);
    try expect(output_tully1D_1.pop.at(0)                   ,  0.57            );
    try expect(output_tully1D_1.pop.at(1)                   ,  0.43            );

    try expect(output_tully1D_2.r.at(0)                     , 21.24034857366279);
    try expect(output_tully1D_2.p.at(0)                     , 20.46292853412406);
    try expect(output_tully1D_2.Ekin + output_tully1D_2.Epot,  0.14951227243131);
    try expect(output_tully1D_2.pop.at(0)                   ,  0.12            );
    try expect(output_tully1D_2.pop.at(1)                   ,  0.88            );

    try expect(output_tully1D_3.r.at(0)                     ,  3.96609063054275);
    try expect(output_tully1D_3.p.at(0)                     ,  1.90698142764801);
    try expect(output_tully1D_3.Ekin + output_tully1D_3.Epot,  0.10014861728192);
    try expect(output_tully1D_3.pop.at(0)                   ,  0.59            );
    try expect(output_tully1D_3.pop.at(1)                   ,  0.41            );
}

test "classical_dynamics_nonadiabatic_lzsh" {
    const opt_tully1D_1 = input.ClassicalDynamicsOptions(f64){
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
        .potential = "tully1D_1",
        .time_step = 1,
        .trajectories = 100,
    };

    var opt_tully1D_2 = opt_tully1D_1; opt_tully1D_2.potential = "tully1D_2";
    var opt_tully1D_3 = opt_tully1D_1; opt_tully1D_3.potential = "tully1D_3";

    const output_tully1D_1 = try classical_dynamics.run(f64, opt_tully1D_1, false, allocator); defer output_tully1D_1.deinit();
    const output_tully1D_2 = try classical_dynamics.run(f64, opt_tully1D_2, false, allocator); defer output_tully1D_2.deinit();
    const output_tully1D_3 = try classical_dynamics.run(f64, opt_tully1D_3, false, allocator); defer output_tully1D_3.deinit();

    try expect(output_tully1D_1.r.at(0)                     ,  21.03442355730371);
    try expect(output_tully1D_1.p.at(0)                     ,  21.11129461625254);
    try expect(output_tully1D_1.Ekin + output_tully1D_1.Epot,   0.10955556586487);
    try expect(output_tully1D_1.pop.at(0)                   ,   0.61            );
    try expect(output_tully1D_1.pop.at(1)                   ,   0.39            );

    try expect(output_tully1D_2.r.at(0)                     ,  21.35503334318419);
    try expect(output_tully1D_2.p.at(0)                     ,  20.66406116097518);
    try expect(output_tully1D_2.Ekin + output_tully1D_2.Epot,   0.14955544536330);
    try expect(output_tully1D_2.pop.at(0)                   ,   0.16            );
    try expect(output_tully1D_2.pop.at(1)                   ,   0.84            );

    try expect(output_tully1D_3.r.at(0)                     , -16.75966919098743);
    try expect(output_tully1D_3.p.at(0)                     , -19.93694370559839);
    try expect(output_tully1D_3.Ekin + output_tully1D_3.Epot,   0.10015338906946);
    try expect(output_tully1D_3.pop.at(0)                   ,   0.00            );
    try expect(output_tully1D_3.pop.at(1)                   ,   1.00            );
}

test "classical_dynamics_nonadiabatic_mash" {
    const opt_tully1D_1 = input.ClassicalDynamicsOptions(f64){
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
        .potential = "tully1D_1",
        .time_step = 1,
        .trajectories = 100,
    };

    var opt_tully1D_2 = opt_tully1D_1; opt_tully1D_2.potential = "tully1D_2";
    var opt_tully1D_3 = opt_tully1D_1; opt_tully1D_3.potential = "tully1D_3";

    const output_tully1D_1 = try classical_dynamics.run(f64, opt_tully1D_1, false, allocator); defer output_tully1D_1.deinit();
    const output_tully1D_2 = try classical_dynamics.run(f64, opt_tully1D_2, false, allocator); defer output_tully1D_2.deinit();
    const output_tully1D_3 = try classical_dynamics.run(f64, opt_tully1D_3, false, allocator); defer output_tully1D_3.deinit();

    try expect(output_tully1D_1.r.at(0)                     , 21.00714185386426);
    try expect(output_tully1D_1.p.at(0)                     , 21.05977392710928);
    try expect(output_tully1D_1.Ekin + output_tully1D_1.Epot,  0.10954664424545);
    try expect(output_tully1D_1.pop.at(0)                   ,  0.59            );
    try expect(output_tully1D_1.pop.at(1)                   ,  0.41            );

    try expect(output_tully1D_2.r.at(0)                     , 21.10923477506220);
    try expect(output_tully1D_2.p.at(0)                     , 20.37931004595746);
    try expect(output_tully1D_2.Ekin + output_tully1D_2.Epot,  0.14951199660272);
    try expect(output_tully1D_2.pop.at(0)                   ,  0.10            );
    try expect(output_tully1D_2.pop.at(1)                   ,  0.90            );

    try expect(output_tully1D_3.r.at(0)                     , -1.14318287393846);
    try expect(output_tully1D_3.p.at(0)                     , -3.60808111892009);
    try expect(output_tully1D_3.Ekin + output_tully1D_3.Epot,  0.10015092475260);
    try expect(output_tully1D_3.pop.at(0)                   ,  0.68            );
    try expect(output_tully1D_3.pop.at(1)                   ,  0.32            );
}
