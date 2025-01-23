const std = @import("std"); const builtin = @import("builtin");

const cdn = @import("acorn").cdn;
const mpt = @import("acorn").mpt;
const qdn = @import("acorn").qdn;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    try std.io.getStdOut().writer().print("ZIG VERSION: {}\n", .{builtin.zig_version});

    const potentials: [8][]const u8 = [_][]const u8{
        "tully1D_1", "tully1D_2", "tully1D_3", "doubleState1D_1", "doubleState1D_2", "tripleState1D_1", "tripleState1D_2", "tripleState1D_3"
    };

    var opt_exact = qdn.QuantumDynamicsOptions(f64){
        .log_intervals = .{
            .iteration = 50
        },
        .grid = .{
            .limits = &[_]f64{-48, 144}, .points = 4096
        },
        .initial_conditions = .{
            .position = &[_]f64{-15}, .momentum = &[_]f64{10}, .gamma = 2, .state = 0, .mass = 2000
        },

        .adiabatic = true,
        .iterations = 500,
        .mode = &[_]u32{0, 1},
        .time_step = 10,

        .potential = "",
    };

    var opt_classic = cdn.ClassicalDynamicsOptions(f64){
        .log_intervals = .{
            .trajectory = 100, .iteration = 500
        },

        .initial_conditions = .{
            .position_mean = &[_]f64{-15.0},
            .position_std  = &[_]f64{0.5},
            .momentum_mean = &[_]f64{10.0},
            .momentum_std  = &[_]f64{1.0},
            .state = 0, .mass = 2000
        },

        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 5000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 1000,

        .potential = "",
    };

    var opt_diabatic_potential = mpt.ModelPotentialOptions(f64){
        .adiabatic = false,
        .limits = &[_]f64{-16, 16},
        .points = 1024,
    };

    var opt_adiabatic_potential = mpt.ModelPotentialOptions(f64){
        .adiabatic = true,
        .limits = &[_]f64{-16, 16},
        .points = 1024,
    };

    var p: u32 = 10; while (p <= 50) : (p += 10) for (potentials) |potential| {

             if (std.mem.eql(u8, potential, "tully1D_1"      )) {opt_exact.initial_conditions.state = 1; opt_classic.initial_conditions.state = 1;}
        else if (std.mem.eql(u8, potential, "tully1D_2"      )) {opt_exact.initial_conditions.state = 1; opt_classic.initial_conditions.state = 1;}
        else if (std.mem.eql(u8, potential, "tully1D_3"      )) {opt_exact.initial_conditions.state = 0; opt_classic.initial_conditions.state = 1;}
        else if (std.mem.eql(u8, potential, "doubleState1D_1")) {opt_exact.initial_conditions.state = 1; opt_classic.initial_conditions.state = 1;}
        else if (std.mem.eql(u8, potential, "doubleState1D_2")) {opt_exact.initial_conditions.state = 1; opt_classic.initial_conditions.state = 1;}
        else if (std.mem.eql(u8, potential, "tripleState1D_1")) {opt_exact.initial_conditions.state = 2; opt_classic.initial_conditions.state = 2;}
        else if (std.mem.eql(u8, potential, "tripleState1D_2")) {opt_exact.initial_conditions.state = 2; opt_classic.initial_conditions.state = 2;}
        else if (std.mem.eql(u8, potential, "tripleState1D_3")) {opt_exact.initial_conditions.state = 1; opt_classic.initial_conditions.state = 2;}

        else return error.StateNotDefined; var opt_fssh = opt_classic; var opt_kfssh = opt_classic; var opt_lzsh = opt_classic; var opt_mash = opt_classic;

        opt_fssh.fewest_switches  = .{.quantum_substep = 10, .decoherence_alpha = null  };
        opt_kfssh.fewest_switches = .{.quantum_substep = 10, .decoherence_alpha = null  };
        opt_lzsh.landau_zener     = .{                                                  };
        opt_mash.spin_mapping     = .{.quantum_substep = 10, .sample_bloch_vector = true};

        opt_exact.initial_conditions.momentum      = &[_]f64{@as(f64, @floatFromInt(p))};
        opt_fssh.initial_conditions.momentum_mean  = &[_]f64{@as(f64, @floatFromInt(p))};
        opt_kfssh.initial_conditions.momentum_mean = &[_]f64{@as(f64, @floatFromInt(p))};
        opt_lzsh.initial_conditions.momentum_mean  = &[_]f64{@as(f64, @floatFromInt(p))};
        opt_mash.initial_conditions.momentum_mean  = &[_]f64{@as(f64, @floatFromInt(p))};

        opt_fssh.time_derivative_coupling  = "numeric";
        opt_kfssh.time_derivative_coupling = "baeckan";
        opt_mash.time_derivative_coupling  = "numeric";

        const  diabatic_potential = try allocator.alloc(u8, 23 + potential.len); defer allocator.free(diabatic_potential);
        const adiabatic_potential = try allocator.alloc(u8, 24 + potential.len); defer allocator.free(adiabatic_potential);

        _ = try std.fmt.bufPrint(diabatic_potential,  "POTENTIAL_DIABATIC_{s}.mat",  .{potential});
        _ = try std.fmt.bufPrint(adiabatic_potential, "POTENTIAL_ADIABATIC_{s}.mat", .{potential});

        opt_diabatic_potential.output    = diabatic_potential; opt_adiabatic_potential.output    = adiabatic_potential;
        opt_diabatic_potential.potential =          potential; opt_adiabatic_potential.potential =           potential;

        const kinetic_energy_exact   = try allocator.alloc(u8, 30 + potential.len); defer allocator.free(kinetic_energy_exact  );
        const momentum_exact         = try allocator.alloc(u8, 24 + potential.len); defer allocator.free(momentum_exact        );
        const population_exact       = try allocator.alloc(u8, 26 + potential.len); defer allocator.free(population_exact      );
        const position_exact         = try allocator.alloc(u8, 24 + potential.len); defer allocator.free(position_exact        );
        const potential_energy_exact = try allocator.alloc(u8, 32 + potential.len); defer allocator.free(potential_energy_exact);
        const total_energy_exact     = try allocator.alloc(u8, 28 + potential.len); defer allocator.free(total_energy_exact    );

        _ = try std.fmt.bufPrint(kinetic_energy_exact,   "KINETIC_ENERGY_{s}_P={d}_EXACT.mat"  , .{potential, p});
        _ = try std.fmt.bufPrint(momentum_exact,         "MOMENTUM_{s}_P={d}_EXACT.mat"        , .{potential, p});
        _ = try std.fmt.bufPrint(population_exact,       "POPULATION_{s}_P={d}_EXACT.mat"      , .{potential, p});
        _ = try std.fmt.bufPrint(position_exact,         "POSITION_{s}_P={d}_EXACT.mat"        , .{potential, p});
        _ = try std.fmt.bufPrint(potential_energy_exact, "POTENTIAL_ENERGY_{s}_P={d}_EXACT.mat", .{potential, p});
        _ = try std.fmt.bufPrint(total_energy_exact,     "TOTAL_ENERGY_{s}_P={d}_EXACT.mat"    , .{potential, p});

        opt_exact.write.kinetic_energy   =   kinetic_energy_exact;
        opt_exact.write.momentum         =         momentum_exact;
        opt_exact.write.population       =       population_exact;
        opt_exact.write.position         =         position_exact;
        opt_exact.write.potential_energy = potential_energy_exact;
        opt_exact.write.total_energy     =     total_energy_exact;
        opt_exact.potential              =              potential;

        const coefficient_mean_fssh              = try allocator.alloc(u8, 31 + potential.len); defer allocator.free(coefficient_mean_fssh             );
        const kinetic_energy_mean_fssh           = try allocator.alloc(u8, 34 + potential.len); defer allocator.free(kinetic_energy_mean_fssh          );
        const momentum_mean_fssh                 = try allocator.alloc(u8, 28 + potential.len); defer allocator.free(momentum_mean_fssh                );
        const population_mean_fssh               = try allocator.alloc(u8, 30 + potential.len); defer allocator.free(population_mean_fssh              );
        const position_mean_fssh                 = try allocator.alloc(u8, 28 + potential.len); defer allocator.free(position_mean_fssh                );
        const potential_energy_mean_fssh         = try allocator.alloc(u8, 36 + potential.len); defer allocator.free(potential_energy_mean_fssh        );
        const time_derivative_coupling_mean_fssh = try allocator.alloc(u8, 44 + potential.len); defer allocator.free(time_derivative_coupling_mean_fssh);
        const total_energy_mean_fssh             = try allocator.alloc(u8, 32 + potential.len); defer allocator.free(total_energy_mean_fssh            );

        _ = try std.fmt.bufPrint(coefficient_mean_fssh,              "COEFFICIENT_MEAN_{s}_P={d}_FSSH.mat",      .{potential, p});
        _ = try std.fmt.bufPrint(kinetic_energy_mean_fssh,           "KINETIC_ENERGY_MEAN_{s}_P={d}_FSSH.mat",   .{potential, p});
        _ = try std.fmt.bufPrint(momentum_mean_fssh,                 "MOMENTUM_MEAN_{s}_P={d}_FSSH.mat",         .{potential, p});
        _ = try std.fmt.bufPrint(population_mean_fssh,               "POPULATION_MEAN_{s}_P={d}_FSSH.mat",       .{potential, p});
        _ = try std.fmt.bufPrint(position_mean_fssh,                 "POSITION_MEAN_{s}_P={d}_FSSH.mat",         .{potential, p});
        _ = try std.fmt.bufPrint(potential_energy_mean_fssh,         "POTENTIAL_ENERGY_MEAN_{s}_P={d}_FSSH.mat", .{potential, p});
        _ = try std.fmt.bufPrint(time_derivative_coupling_mean_fssh, "TDC_MEAN_{s}_P={d}_FSSH.mat",              .{potential, p});
        _ = try std.fmt.bufPrint(total_energy_mean_fssh,             "TOTAL_ENERGY_MEAN_{s}_P={d}_FSSH.mat",     .{potential, p});

        opt_fssh.write.fssh_coefficient_mean         =              coefficient_mean_fssh;
        opt_fssh.write.kinetic_energy_mean           =           kinetic_energy_mean_fssh;
        opt_fssh.write.momentum_mean                 =                 momentum_mean_fssh;
        opt_fssh.write.population_mean               =               population_mean_fssh;
        opt_fssh.write.position_mean                 =                 position_mean_fssh;
        opt_fssh.write.potential_energy_mean         =         potential_energy_mean_fssh;
        opt_fssh.write.time_derivative_coupling_mean = time_derivative_coupling_mean_fssh;
        opt_fssh.write.total_energy_mean             =             total_energy_mean_fssh;
        opt_fssh.potential                           =                          potential;

        const coefficient_mean_kfssh              = try allocator.alloc(u8, 32 + potential.len); defer allocator.free(coefficient_mean_kfssh             );
        const kinetic_energy_mean_kfssh           = try allocator.alloc(u8, 35 + potential.len); defer allocator.free(kinetic_energy_mean_kfssh          );
        const momentum_mean_kfssh                 = try allocator.alloc(u8, 29 + potential.len); defer allocator.free(momentum_mean_kfssh                );
        const population_mean_kfssh               = try allocator.alloc(u8, 31 + potential.len); defer allocator.free(population_mean_kfssh              );
        const position_mean_kfssh                 = try allocator.alloc(u8, 29 + potential.len); defer allocator.free(position_mean_kfssh                );
        const potential_energy_mean_kfssh         = try allocator.alloc(u8, 37 + potential.len); defer allocator.free(potential_energy_mean_kfssh        );
        const time_derivative_coupling_mean_kfssh = try allocator.alloc(u8, 45 + potential.len); defer allocator.free(time_derivative_coupling_mean_kfssh);
        const total_energy_mean_kfssh             = try allocator.alloc(u8, 33 + potential.len); defer allocator.free(total_energy_mean_kfssh            );

        _ = try std.fmt.bufPrint(coefficient_mean_kfssh,              "COEFFICIENT_MEAN_{s}_P={d}_KFSSH.mat",      .{potential, p});
        _ = try std.fmt.bufPrint(kinetic_energy_mean_kfssh,           "KINETIC_ENERGY_MEAN_{s}_P={d}_KFSSH.mat",   .{potential, p});
        _ = try std.fmt.bufPrint(momentum_mean_kfssh,                 "MOMENTUM_MEAN_{s}_P={d}_KFSSH.mat",         .{potential, p});
        _ = try std.fmt.bufPrint(population_mean_kfssh,               "POPULATION_MEAN_{s}_P={d}_KFSSH.mat",       .{potential, p});
        _ = try std.fmt.bufPrint(position_mean_kfssh,                 "POSITION_MEAN_{s}_P={d}_KFSSH.mat",         .{potential, p});
        _ = try std.fmt.bufPrint(potential_energy_mean_kfssh,         "POTENTIAL_ENERGY_MEAN_{s}_P={d}_KFSSH.mat", .{potential, p});
        _ = try std.fmt.bufPrint(time_derivative_coupling_mean_kfssh, "TDC_MEAN_{s}_P={d}_KFSSH.mat",              .{potential, p});
        _ = try std.fmt.bufPrint(total_energy_mean_kfssh,             "TOTAL_ENERGY_MEAN_{s}_P={d}_KFSSH.mat",     .{potential, p});

        opt_kfssh.write.fssh_coefficient_mean         =              coefficient_mean_kfssh;
        opt_kfssh.write.kinetic_energy_mean           =           kinetic_energy_mean_kfssh;
        opt_kfssh.write.momentum_mean                 =                 momentum_mean_kfssh;
        opt_kfssh.write.population_mean               =               population_mean_kfssh;
        opt_kfssh.write.position_mean                 =                 position_mean_kfssh;
        opt_kfssh.write.potential_energy_mean         =         potential_energy_mean_kfssh;
        opt_kfssh.write.time_derivative_coupling_mean = time_derivative_coupling_mean_kfssh;
        opt_kfssh.write.total_energy_mean             =             total_energy_mean_kfssh;
        opt_kfssh.potential                           =                           potential;

        const kinetic_energy_mean_lzsh   = try allocator.alloc(u8, 34 + potential.len); defer allocator.free(kinetic_energy_mean_lzsh  );
        const momentum_mean_lzsh         = try allocator.alloc(u8, 28 + potential.len); defer allocator.free(momentum_mean_lzsh        );
        const population_mean_lzsh       = try allocator.alloc(u8, 30 + potential.len); defer allocator.free(population_mean_lzsh      );
        const position_mean_lzsh         = try allocator.alloc(u8, 28 + potential.len); defer allocator.free(position_mean_lzsh        );
        const potential_energy_mean_lzsh = try allocator.alloc(u8, 36 + potential.len); defer allocator.free(potential_energy_mean_lzsh);
        const total_energy_mean_lzsh     = try allocator.alloc(u8, 32 + potential.len); defer allocator.free(total_energy_mean_lzsh    );

        _ = try std.fmt.bufPrint(kinetic_energy_mean_lzsh,   "KINETIC_ENERGY_MEAN_{s}_P={d}_LZSH.mat",   .{potential, p});
        _ = try std.fmt.bufPrint(momentum_mean_lzsh,         "MOMENTUM_MEAN_{s}_P={d}_LZSH.mat",         .{potential, p});
        _ = try std.fmt.bufPrint(population_mean_lzsh,       "POPULATION_MEAN_{s}_P={d}_LZSH.mat",       .{potential, p});
        _ = try std.fmt.bufPrint(position_mean_lzsh,         "POSITION_MEAN_{s}_P={d}_LZSH.mat",         .{potential, p});
        _ = try std.fmt.bufPrint(potential_energy_mean_lzsh, "POTENTIAL_ENERGY_MEAN_{s}_P={d}_LZSH.mat", .{potential, p});
        _ = try std.fmt.bufPrint(total_energy_mean_lzsh,     "TOTAL_ENERGY_MEAN_{s}_P={d}_LZSH.mat",     .{potential, p});

        opt_lzsh.write.kinetic_energy_mean           =           kinetic_energy_mean_lzsh;
        opt_lzsh.write.momentum_mean                 =                 momentum_mean_lzsh;
        opt_lzsh.write.population_mean               =               population_mean_lzsh;
        opt_lzsh.write.position_mean                 =                 position_mean_lzsh;
        opt_lzsh.write.potential_energy_mean         =         potential_energy_mean_lzsh;
        opt_lzsh.write.total_energy_mean             =             total_energy_mean_lzsh;
        opt_lzsh.potential                           =                          potential;

        const kinetic_energy_mean_mash           = try allocator.alloc(u8, 34 + potential.len); defer allocator.free(kinetic_energy_mean_mash          );
        const momentum_mean_mash                 = try allocator.alloc(u8, 28 + potential.len); defer allocator.free(momentum_mean_mash                );
        const population_mean_mash               = try allocator.alloc(u8, 30 + potential.len); defer allocator.free(population_mean_mash              );
        const position_mean_mash                 = try allocator.alloc(u8, 28 + potential.len); defer allocator.free(position_mean_mash                );
        const potential_energy_mean_mash         = try allocator.alloc(u8, 36 + potential.len); defer allocator.free(potential_energy_mean_mash        );
        const time_derivative_coupling_mean_mash = try allocator.alloc(u8, 44 + potential.len); defer allocator.free(time_derivative_coupling_mean_mash);
        const total_energy_mean_mash             = try allocator.alloc(u8, 32 + potential.len); defer allocator.free(total_energy_mean_mash            );

        _ = try std.fmt.bufPrint(kinetic_energy_mean_mash,           "KINETIC_ENERGY_MEAN_{s}_P={d}_MASH.mat",   .{potential, p});
        _ = try std.fmt.bufPrint(momentum_mean_mash,                 "MOMENTUM_MEAN_{s}_P={d}_MASH.mat",         .{potential, p});
        _ = try std.fmt.bufPrint(population_mean_mash,               "POPULATION_MEAN_{s}_P={d}_MASH.mat",       .{potential, p});
        _ = try std.fmt.bufPrint(position_mean_mash,                 "POSITION_MEAN_{s}_P={d}_MASH.mat",         .{potential, p});
        _ = try std.fmt.bufPrint(potential_energy_mean_mash,         "POTENTIAL_ENERGY_MEAN_{s}_P={d}_MASH.mat", .{potential, p});
        _ = try std.fmt.bufPrint(time_derivative_coupling_mean_mash, "TDC_MEAN_{s}_P={d}_MASH.mat",              .{potential, p});
        _ = try std.fmt.bufPrint(total_energy_mean_mash,             "TOTAL_ENERGY_MEAN_{s}_P={d}_MASH.mat",     .{potential, p});

        opt_mash.write.kinetic_energy_mean           =           kinetic_energy_mean_mash;
        opt_mash.write.momentum_mean                 =                 momentum_mean_mash;
        opt_mash.write.population_mean               =               population_mean_mash;
        opt_mash.write.position_mean                 =                 position_mean_mash;
        opt_mash.write.potential_energy_mean         =         potential_energy_mean_mash;
        opt_mash.write.time_derivative_coupling_mean = time_derivative_coupling_mean_mash;
        opt_mash.write.total_energy_mean             =             total_energy_mean_mash;
        opt_mash.potential                           =                          potential;

        try mpt.write(f64, opt_diabatic_potential, allocator); try mpt.write(f64, opt_adiabatic_potential, allocator);

        const output_exact = try qdn.run(f64, opt_exact, true, allocator); defer output_exact.deinit();
        const output_fssh  = try cdn.run(f64, opt_fssh , true, allocator); defer  output_fssh.deinit();
        const output_kfssh = try cdn.run(f64, opt_kfssh, true, allocator); defer output_kfssh.deinit();
        const output_lzsh  = try cdn.run(f64, opt_lzsh , true, allocator); defer  output_lzsh.deinit();

        if (!std.mem.eql(u8, potential[0..6], "triple")) {
            const output_mash  = try cdn.run(f64, opt_mash , true, allocator); defer output_mash.deinit();
        }
    };
}
