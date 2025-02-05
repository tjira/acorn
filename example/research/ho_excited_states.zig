const std = @import("std"); const builtin = @import("builtin");

const quantum_dynamics = @import("acorn").quantum_dynamics;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    try std.io.getStdOut().writer().print("ZIG VERSION: {}\n", .{builtin.zig_version});

    const opt_imaginary_1D_1 = quantum_dynamics.QuantumDynamicsOptions(f64){
        .log_intervals = .{
            .iteration = 200
        },

        .grid = .{
            .limits = &[_]f64{-8, 8}, .points = 64
        },

        .initial_conditions = .{
            .position = &[_]f64{1}, .momentum = &[_]f64{0}, .gamma = 2, .state = 0, .mass = 1
        },

        .write = .{
            // .wavefunction = "WAVEFUNCTION_IMAG_1D_1.mat",
        },

        .adiabatic = false,
        .iterations = 1000,
        .mode = &[_]u32{1, 0},
        .potential = "harmonic1D_1",
        .time_step = 0.1,
    };

    const opt_imaginary_2D_1 = quantum_dynamics.QuantumDynamicsOptions(f64){
        .log_intervals = .{
            .iteration = 200
        },

        .grid = .{
            .limits = &[_]f64{-8, 8}, .points = 64
        },

        .initial_conditions = .{
            .position = &[_]f64{1, 1}, .momentum = &[_]f64{0, 0}, .gamma = 2, .state = 0, .mass = 1
        },

        .write = .{
            // .wavefunction = "WAVEFUNCTION_IMAG_2D_1.mat",
        },

        .adiabatic = false,
        .iterations = 1000,
        .mode = &[_]u32{1, 0},
        .potential = "harmonic2D_1",
        .time_step = 0.1,
    };

    const opt_imaginary_3D_1 = quantum_dynamics.QuantumDynamicsOptions(f64){
        .log_intervals = .{
            .iteration = 200
        },

        .grid = .{
            .limits = &[_]f64{-8, 8}, .points = 64
        },

        .initial_conditions = .{
            .position = &[_]f64{1, 1, 1}, .momentum = &[_]f64{0, 0, 0}, .gamma = 2, .state = 0, .mass = 1
        },

        .write = .{
            // .wavefunction = "WAVEFUNCTION_IMAG_3D_1.mat",
        },

        .adiabatic = false,
        .iterations = 1000,
        .mode = &[_]u32{1, 0},
        .potential = "harmonic3D_1",
        .time_step = 0.1,
    };

    var opt_imaginary_1D_2 = opt_imaginary_1D_1; opt_imaginary_1D_2.mode = &[_]u32{2, 0}; //opt_imaginary_1D_2.write.wavefunction = "WAVEFUNCTION_IMAG_1D_2.mat";
    var opt_imaginary_1D_3 = opt_imaginary_1D_1; opt_imaginary_1D_3.mode = &[_]u32{3, 0}; //opt_imaginary_1D_3.write.wavefunction = "WAVEFUNCTION_IMAG_1D_3.mat";
    var opt_imaginary_1D_4 = opt_imaginary_1D_1; opt_imaginary_1D_4.mode = &[_]u32{4, 0}; //opt_imaginary_1D_4.write.wavefunction = "WAVEFUNCTION_IMAG_1D_4.mat";
    var opt_imaginary_1D_5 = opt_imaginary_1D_1; opt_imaginary_1D_5.mode = &[_]u32{5, 0}; //opt_imaginary_1D_5.write.wavefunction = "WAVEFUNCTION_IMAG_1D_5.mat";
    var opt_imaginary_1D_6 = opt_imaginary_1D_1; opt_imaginary_1D_6.mode = &[_]u32{6, 0}; //opt_imaginary_1D_6.write.wavefunction = "WAVEFUNCTION_IMAG_1D_6.mat";
    var opt_imaginary_1D_7 = opt_imaginary_1D_1; opt_imaginary_1D_7.mode = &[_]u32{7, 0}; //opt_imaginary_1D_7.write.wavefunction = "WAVEFUNCTION_IMAG_1D_7.mat";
    var opt_imaginary_1D_8 = opt_imaginary_1D_1; opt_imaginary_1D_8.mode = &[_]u32{8, 0}; //opt_imaginary_1D_8.write.wavefunction = "WAVEFUNCTION_IMAG_1D_8.mat";
    var opt_imaginary_1D_9 = opt_imaginary_1D_1; opt_imaginary_1D_9.mode = &[_]u32{9, 0}; //opt_imaginary_1D_9.write.wavefunction = "WAVEFUNCTION_IMAG_1D_9.mat";

    var opt_imaginary_2D_2 = opt_imaginary_2D_1; opt_imaginary_2D_2.mode = &[_]u32{2, 0}; //opt_imaginary_2D_2.write.wavefunction = "WAVEFUNCTION_IMAG_2D_2.mat";
    var opt_imaginary_2D_3 = opt_imaginary_2D_1; opt_imaginary_2D_3.mode = &[_]u32{3, 0}; //opt_imaginary_2D_3.write.wavefunction = "WAVEFUNCTION_IMAG_2D_3.mat";
    var opt_imaginary_2D_4 = opt_imaginary_2D_1; opt_imaginary_2D_4.mode = &[_]u32{4, 0}; //opt_imaginary_2D_4.write.wavefunction = "WAVEFUNCTION_IMAG_2D_4.mat";
    var opt_imaginary_2D_5 = opt_imaginary_2D_1; opt_imaginary_2D_5.mode = &[_]u32{5, 0}; //opt_imaginary_2D_5.write.wavefunction = "WAVEFUNCTION_IMAG_2D_5.mat";
    var opt_imaginary_2D_6 = opt_imaginary_2D_1; opt_imaginary_2D_6.mode = &[_]u32{6, 0}; //opt_imaginary_2D_6.write.wavefunction = "WAVEFUNCTION_IMAG_2D_6.mat";
    var opt_imaginary_2D_7 = opt_imaginary_2D_1; opt_imaginary_2D_7.mode = &[_]u32{7, 0}; //opt_imaginary_2D_7.write.wavefunction = "WAVEFUNCTION_IMAG_2D_7.mat";
    var opt_imaginary_2D_8 = opt_imaginary_2D_1; opt_imaginary_2D_8.mode = &[_]u32{8, 0}; //opt_imaginary_2D_8.write.wavefunction = "WAVEFUNCTION_IMAG_2D_8.mat";
    var opt_imaginary_2D_9 = opt_imaginary_2D_1; opt_imaginary_2D_9.mode = &[_]u32{9, 0}; //opt_imaginary_2D_9.write.wavefunction = "WAVEFUNCTION_IMAG_2D_9.mat";

    var opt_real_1D = opt_imaginary_1D_1; opt_real_1D.mode = &[_]u32{0, 1};
    var opt_real_2D = opt_imaginary_2D_1; opt_real_2D.mode = &[_]u32{0, 1};
    var opt_real_3D = opt_imaginary_3D_1; opt_real_3D.mode = &[_]u32{0, 1};

    opt_real_1D.write.spectrum = "SPECTRUM_1D.mat";
    opt_real_2D.write.spectrum = "SPECTRUM_2D.mat";
    opt_real_3D.write.spectrum = "SPECTRUM_3D.mat";

    opt_real_1D.write.autocorrelation_function = "ACF_1D.mat";
    opt_real_2D.write.autocorrelation_function = "ACF_2D.mat";
    opt_real_3D.write.autocorrelation_function = "ACF_3D.mat";

    // opt_real_1D.write.wavefunction = "WAVEFUNCTION_1D_REAL.mat";
    // opt_real_2D.write.wavefunction = "WAVEFUNCTION_2D_REAL.mat";
    // opt_real_3D.write.wavefunction = "WAVEFUNCTION_3D_REAL.mat";

    const output_imaginary_1D_1 = try quantum_dynamics.run(f64, opt_imaginary_1D_1, true, allocator); defer output_imaginary_1D_1.deinit();
    const output_imaginary_1D_2 = try quantum_dynamics.run(f64, opt_imaginary_1D_2, true, allocator); defer output_imaginary_1D_2.deinit();
    const output_imaginary_1D_3 = try quantum_dynamics.run(f64, opt_imaginary_1D_3, true, allocator); defer output_imaginary_1D_3.deinit();
    const output_imaginary_1D_4 = try quantum_dynamics.run(f64, opt_imaginary_1D_4, true, allocator); defer output_imaginary_1D_4.deinit();
    const output_imaginary_1D_5 = try quantum_dynamics.run(f64, opt_imaginary_1D_5, true, allocator); defer output_imaginary_1D_5.deinit();
    const output_imaginary_1D_6 = try quantum_dynamics.run(f64, opt_imaginary_1D_6, true, allocator); defer output_imaginary_1D_6.deinit();
    const output_imaginary_1D_7 = try quantum_dynamics.run(f64, opt_imaginary_1D_7, true, allocator); defer output_imaginary_1D_7.deinit();
    const output_imaginary_1D_8 = try quantum_dynamics.run(f64, opt_imaginary_1D_8, true, allocator); defer output_imaginary_1D_8.deinit();
    const output_imaginary_1D_9 = try quantum_dynamics.run(f64, opt_imaginary_1D_9, true, allocator); defer output_imaginary_1D_9.deinit();
    const output_real_1D        = try quantum_dynamics.run(f64, opt_real_1D,        true, allocator); defer        output_real_1D.deinit();

    const output_imaginary_2D_1 = try quantum_dynamics.run(f64, opt_imaginary_2D_1, true, allocator); defer output_imaginary_2D_1.deinit();
    const output_imaginary_2D_2 = try quantum_dynamics.run(f64, opt_imaginary_2D_2, true, allocator); defer output_imaginary_2D_2.deinit();
    const output_imaginary_2D_3 = try quantum_dynamics.run(f64, opt_imaginary_2D_3, true, allocator); defer output_imaginary_2D_3.deinit();
    const output_imaginary_2D_4 = try quantum_dynamics.run(f64, opt_imaginary_2D_4, true, allocator); defer output_imaginary_2D_4.deinit();
    const output_imaginary_2D_5 = try quantum_dynamics.run(f64, opt_imaginary_2D_5, true, allocator); defer output_imaginary_2D_5.deinit();
    const output_imaginary_2D_6 = try quantum_dynamics.run(f64, opt_imaginary_2D_6, true, allocator); defer output_imaginary_2D_6.deinit();
    const output_imaginary_2D_7 = try quantum_dynamics.run(f64, opt_imaginary_2D_7, true, allocator); defer output_imaginary_2D_7.deinit();
    const output_imaginary_2D_8 = try quantum_dynamics.run(f64, opt_imaginary_2D_8, true, allocator); defer output_imaginary_2D_8.deinit();
    const output_imaginary_2D_9 = try quantum_dynamics.run(f64, opt_imaginary_2D_9, true, allocator); defer output_imaginary_2D_9.deinit();
    const output_real_2D        = try quantum_dynamics.run(f64, opt_real_2D,        true, allocator); defer        output_real_2D.deinit();

    const output_real_3D = try quantum_dynamics.run(f64, opt_real_3D, true, allocator); defer output_real_3D.deinit();
}
