const std = @import("std"); const builtin = @import("builtin");

const quantum_dynamics = @import("acorn").quantum_dynamics;
const model_potential  = @import("acorn").model_potential ;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    try std.io.getStdOut().writer().print("ZIG VERSION: {}\n", .{builtin.zig_version});

    const opt_diabatic_tully1D_1 = quantum_dynamics.QuantumDynamicsOptions(f64){
        .log_intervals = .{
            .iteration = 50
        },

        .grid = .{
            .limits = &[_]f64{-24, 36}, .points = 4096
        },

        .initial_conditions = .{
            .position = &[_]f64{-15}, .momentum = &[_]f64{35}, .gamma = 2, .state = 1, .mass = 2000
        },

        .write = .{
            .population = "POPULATION_DIABATIC_tully1D_1.mat",
            .wavefunction = "WAVEFUNCTION_DIABATIC_tully1D_1.mat"
        },

        .adiabatic = false,
        .iterations = 200,
        .mode = &[_]u32{0, 1},
        .potential = "tully1D_1",
        .time_step = 10,
    };

    const opt_diabatic_potential_tully1D_1 = model_potential.ModelPotentialOptions(f64){
        .adiabatic = false,
        .limits = &[_]f64{-24, 36},
        .output = "POTENTIAL_DIABATIC_tully1D_1.mat",
        .points = 4096,
        .potential = "tully1D_1"
    };

    var opt_adiabatic_tully1D_1 = opt_diabatic_tully1D_1; opt_adiabatic_tully1D_1.adiabatic = true;
    opt_adiabatic_tully1D_1.write = .{
        .population = "POPULATION_ADIABATIC_tully1D_1.mat",
        .wavefunction = "WAVEFUNCTION_ADIABATIC_tully1D_1.mat"
    };

    var opt_adiabatic_potential_tully1D_1 = opt_diabatic_potential_tully1D_1;
    opt_adiabatic_potential_tully1D_1.adiabatic = true;
    opt_adiabatic_potential_tully1D_1.output = "POTENTIAL_ADIABATIC_tully1D_1.mat";

    const output_diabatic_tully1D_1  = try quantum_dynamics.run(f64, opt_diabatic_tully1D_1,  true, allocator); defer  output_diabatic_tully1D_1.deinit();
    const output_adiabatic_tully1D_1 = try quantum_dynamics.run(f64, opt_adiabatic_tully1D_1, true, allocator); defer output_adiabatic_tully1D_1.deinit();

    try model_potential.write(f64, opt_diabatic_potential_tully1D_1,  allocator);
    try model_potential.write(f64, opt_adiabatic_potential_tully1D_1, allocator);
}
