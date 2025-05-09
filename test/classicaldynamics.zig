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

    try expect(output_harmonic1D_1.r.at(0),                         -0.90963222375214);
    try expect(output_harmonic1D_1.p.at(0),                          0.54149007140784);
    try expect(output_harmonic1D_1.Ekin + output_harmonic1D_1.Epot,  1.22409230923698);
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

    try expect(output_tully1D_1.r.at(0)                     , 20.98607919174036);
    try expect(output_tully1D_1.p.at(0)                     , 21.04378573691536);
    try expect(output_tully1D_1.Ekin + output_tully1D_1.Epot,  0.10955046482012);
    try expect(output_tully1D_1.pop.at(0)                   ,  0.58            );
    try expect(output_tully1D_1.pop.at(1)                   ,  0.42            );
    try expect(output_tully1D_2.r.at(0)                     , 22.52319843672515);
    try expect(output_tully1D_2.p.at(0)                     , 21.93331647288330);
    try expect(output_tully1D_2.Ekin + output_tully1D_2.Epot,  0.14964269361075);
    try expect(output_tully1D_2.pop.at(0)                   ,  0.44            );
    try expect(output_tully1D_2.pop.at(1)                   ,  0.56            );
    try expect(output_tully1D_3.r.at(0)                     ,  5.53685701300704);
    try expect(output_tully1D_3.p.at(0)                     ,  3.54174478703978);
    try expect(output_tully1D_3.Ekin + output_tully1D_3.Epot,  0.10095896173608);
    try expect(output_tully1D_3.pop.at(0)                   ,  0.57            );
    try expect(output_tully1D_3.pop.at(1)                   ,  0.43            );
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

    try expect(output_tully1D_1.r.at(0)                     ,  21.03423896303616);
    try expect(output_tully1D_1.p.at(0)                     ,  21.11141112677523);
    try expect(output_tully1D_1.Ekin + output_tully1D_1.Epot,   0.10955678747866);
    try expect(output_tully1D_1.pop.at(0)                   ,   0.61            );
    try expect(output_tully1D_1.pop.at(1)                   ,   0.39            );

    try expect(output_tully1D_2.r.at(0)                     ,  21.36247947689098);
    try expect(output_tully1D_2.p.at(0)                     ,  20.67185596545007);
    try expect(output_tully1D_2.Ekin + output_tully1D_2.Epot,   0.14963564613228);
    try expect(output_tully1D_2.pop.at(0)                   ,   0.16            );
    try expect(output_tully1D_2.pop.at(1)                   ,   0.84            );

    try expect(output_tully1D_3.r.at(0)                     , -16.84693961792927);
    try expect(output_tully1D_3.p.at(0)                     , -20.05554167860687);
    try expect(output_tully1D_3.Ekin + output_tully1D_3.Epot,   0.10134219509159);
    try expect(output_tully1D_3.pop.at(0)                   ,   0.0             );
    try expect(output_tully1D_3.pop.at(1)                   ,   1.0             );
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

    try expect(output_tully1D_1.r.at(0)                     , 21.00695512171220);
    try expect(output_tully1D_1.p.at(0)                     , 21.05989704416699);
    try expect(output_tully1D_1.Ekin + output_tully1D_1.Epot,  0.10954793322135);
    try expect(output_tully1D_1.pop.at(0)                   ,  0.59            );
    try expect(output_tully1D_1.pop.at(1)                   ,  0.41            );
    try expect(output_tully1D_2.r.at(0)                     , 22.10044168360206);
    try expect(output_tully1D_2.p.at(0)                     , 21.45661238350094);
    try expect(output_tully1D_2.Ekin + output_tully1D_2.Epot,  0.14967491941739);
    try expect(output_tully1D_2.pop.at(0)                   ,  0.33            );
    try expect(output_tully1D_2.pop.at(1)                   ,  0.67            );
    try expect(output_tully1D_3.r.at(0)                     , -1.19219562810638);
    try expect(output_tully1D_3.p.at(0)                     , -3.66930755012664);
    try expect(output_tully1D_3.Ekin + output_tully1D_3.Epot,  0.10108100457442);
    try expect(output_tully1D_3.pop.at(0)                   ,  0.55            );
    try expect(output_tully1D_3.pop.at(1)                   ,  0.45            );
}
