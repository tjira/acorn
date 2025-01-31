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
            .state = &[_]f64{1}, .mass = &[_]f64{1}
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

    try std.testing.expect(@abs(output_tully1D_1.r.at(0)                      - 20.95696749749547) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_1.p.at(0)                      - 21.01999888560924) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_1.Ekin + output_tully1D_1.Epot -  0.10954931488450) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_1.pop.at(0)                    -  0.57            ) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_1.pop.at(1)                    -  0.43            ) < 1e-14);

    try std.testing.expect(@abs(output_tully1D_2.r.at(0)                      - 21.23011710939573) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_2.p.at(0)                      - 20.46292853412406) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_2.Ekin + output_tully1D_2.Epot -  0.14951227243131) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_2.pop.at(0)                    -  0.12            ) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_2.pop.at(1)                    -  0.88            ) < 1e-14);

    try std.testing.expect(@abs(output_tully1D_3.r.at(0)                      - 3.96513713982892) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_3.p.at(0)                      - 1.90698142761924) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_3.Ekin + output_tully1D_3.Epot - 0.10014861728193) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_3.pop.at(0)                    - 0.59            ) < 1e-14);
    try std.testing.expect(@abs(output_tully1D_3.pop.at(1)                    - 0.41            ) < 1e-14);
}
