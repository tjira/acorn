const std = @import("std");

const classical_dynamics = @import("acorn").classical_dynamics;

const allocator = std.testing.allocator;

test "classical_dynamics_adiabatic" {
    const opt_harmonic1D_1 = classical_dynamics.ClassicalDynamicsOptions(f64){
        .initial_conditions = .{
            .position_mean = &[_]f64{1.0},
            .position_std  = &[_]f64{0.5},
            .momentum_mean = &[_]f64{0.0},
            .momentum_std  = &[_]f64{1.0},
            .state = 0, .mass = &[_]f64{1}
        },

        .adiabatic = true,
        .iterations = 1000,
        .potential = "harmonic1D_1",
        .time_step = 0.01,
        .trajectories = 1000,
    };

    const output_harmonic1D_1 = try classical_dynamics.run(f64, opt_harmonic1D_1, false, allocator); defer output_harmonic1D_1.deinit();

    try std.testing.expect(@abs(output_harmonic1D_1.r.at(0)                         + 0.86752945241864) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.p.at(0)                         - 0.50663328133126) < 1e-14);
    try std.testing.expect(@abs(output_harmonic1D_1.Ekin + output_harmonic1D_1.Epot - 1.10695573782275) < 1e-14);
}

test "classical_dynamics_nonadiabatic_fssh" {
    const opt_tully1D_1 = classical_dynamics.ClassicalDynamicsOptions(f64){
        .initial_conditions = .{
            .position_mean = &[_]f64{-15.0},
            .position_std  = &[_]f64{  0.5},
            .momentum_mean = &[_]f64{ 15.0},
            .momentum_std  = &[_]f64{  1.0},
            .state = 0, .mass = &[_]f64{2000}
        },

        .fewest_switches = .{},

        .adiabatic = true,
        .iterations = 5000,
        .potential = "tully1D_1",
        .time_step = 1,
        .trajectories = 100,
    };

    const output_tully1D_1 = try classical_dynamics.run(f64, opt_tully1D_1, false, allocator); defer output_tully1D_1.deinit();

    try std.testing.expect(@abs(output_tully1D_1.r.at(0)                      - 20.83715501229083) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_1.p.at(0)                      - 14.01201053460730) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_1.Ekin + output_tully1D_1.Epot -  0.04596357032711) < 1e-14);

    try std.testing.expect(@abs(output_tully1D_1.pop.at(0) - 0.69) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_1.pop.at(1) - 0.31) < 1e-14);
}
