//! Main file of the program.

const std = @import("std"); const builtin = @import("builtin"); const Complex = std.math.complex.Complex;

const allocator = std.heap.page_allocator; const fsize = 2048;

pub const basis                     = @import("basis.zig"                   );
pub const blas                      = @import("blas.zig"                    );
pub const classical_dynamics        = @import("classicaldynamics.zig"       );
pub const configuration_interaction = @import("configurationinteraction.zig");
pub const constant                  = @import("constant.zig"                );
pub const contracted_gaussian       = @import("contractedgaussian.zig"      );
pub const fftw                      = @import("fftw.zig"                    );
pub const fibonacci                 = @import("fibonacci.zig"               );
pub const hartree_fock              = @import("hartreefock.zig"             );
pub const helper                    = @import("helper.zig"                  );
pub const input                     = @import("input.zig"                   );
pub const integral                  = @import("integral.zig"                );
pub const lapack                    = @import("lapack.zig"                  );
pub const linear_algebra            = @import("linearalgebra.zig"           );
pub const math                      = @import("math.zig"                    );
pub const matrix                    = @import("matrix.zig"                  );
pub const model_potential           = @import("modelpotential.zig"          );
pub const moller_plesset            = @import("mollerplesset.zig"           );
pub const output                    = @import("input.zig"                   );
pub const prime                     = @import("prime.zig"                   );
pub const primitive_gaussian        = @import("primitivegaussian.zig"       );
pub const quantum_dynamics          = @import("quantumdynamics.zig"         );
pub const sort                      = @import("sort.zig"                    );
pub const strided_array             = @import("stridedarray.zig"            );
pub const system                    = @import("system.zig"                  );
pub const tensor                    = @import("tensor.zig"                  );
pub const transform                 = @import("transform.zig"               );
pub const vector                    = @import("vector.zig"                  );
pub const wavefunction              = @import("wavefunction.zig"            );

pub const Matrix = @import("matrix.zig").Matrix;
pub const Tensor = @import("tensor.zig").Tensor;
pub const Vector = @import("vector.zig").Vector;

/// Parse the input JSON file and run the corresponding target.
pub fn parse(filebuf: []const u8) !void {
    const inputjs = try std.json.parseFromSlice(std.json.Value, allocator, filebuf, .{}); defer inputjs.deinit();

    if (inputjs.value.object.contains("classical_dynamics")) {

        const obj = try std.json.parseFromValue(input.ClassicalDynamicsOptions(f64), allocator, inputjs.value.object.get("classical_dynamics").?, .{}); defer obj.deinit();

        _ = try classical_dynamics.run(f64, obj.value, true, allocator);
    }

    if (inputjs.value.object.contains("configuration_interaction")) {

        const obj = try std.json.parseFromValue(input.ConfigurationInteractionOptions(f64), allocator, inputjs.value.object.get("configuration_interaction").?, .{}); defer obj.deinit();

        _ = try configuration_interaction.run(f64, obj.value, true, allocator);
    }

    if (inputjs.value.object.contains("fibonacci")) {

        const obj = try std.json.parseFromValue(input.FibonacciOptions(u128), allocator, inputjs.value.object.get("fibonacci").?, .{}); defer obj.deinit();

        _ = try fibonacci.run(u128, obj.value, true, allocator);
    }

    if (inputjs.value.object.contains("hartree_fock")) {

        const obj = try std.json.parseFromValue(input.HartreeFockOptions(f64), allocator, inputjs.value.object.get("hartree_fock").?, .{}); defer obj.deinit();

        _ = try hartree_fock.run(f64, obj.value, true, allocator);
    }

    if (inputjs.value.object.contains("model_potential")) {

        const obj = try std.json.parseFromValue(input.ModelPotentialOptions(f64), allocator, inputjs.value.object.get("model_potential").?, .{}); defer obj.deinit();

        _ = try model_potential.write(f64, obj.value, allocator);
    }

    if (inputjs.value.object.contains("moller_plesset")) {

        const obj = try std.json.parseFromValue(input.MollerPlessetOptions(f64), allocator, inputjs.value.object.get("moller_plesset").?, .{}); defer obj.deinit();

        _ = try moller_plesset.run(f64, obj.value, true, allocator);
    }

    if (inputjs.value.object.contains("prime")) {

        const obj = try std.json.parseFromValue(input.PrimeOptions(u128), allocator, inputjs.value.object.get("prime").?, .{}); defer obj.deinit();

        _ = try prime.run(u128, obj.value, true, allocator);
    }

    if (inputjs.value.object.contains("quantum_dynamics")) {

        const obj = try std.json.parseFromValue(input.QuantumDynamicsOptions(f64), allocator, inputjs.value.object.get("quantum_dynamics").?, .{}); defer obj.deinit();

        _ = try quantum_dynamics.run(f64, obj.value, true, allocator);
    }

    if (inputjs.value.object.contains("sort")) {

        const obj = try std.json.parseFromValue(input.SortOptions(), allocator, inputjs.value.object.get("sort").?, .{}); defer obj.deinit();

        _ = try sort.run(f64, obj.value, true, allocator);
    }
}

/// The main function of the program.
pub fn main() !void {
    var timer = try std.time.Timer.start(); var argv = try std.process.argsWithAllocator(allocator); defer argv.deinit(); var argc: usize = 0;

    try std.io.getStdOut().writer().print("ZIG VERSION: {}\n", .{builtin.zig_version});

    _ = argv.next(); while (argv.next()) |arg| {

        try std.io.getStdOut().writer().print("\nPROCESSED FILE: {s}\n", .{arg});

        const filebuf = try std.fs.cwd().readFileAlloc(allocator, arg, fsize); defer allocator.free(filebuf);

        try parse(filebuf); argc += 1;
    }

    default: {if (argc == 0) {

        const filebuf = std.fs.cwd().readFileAlloc(allocator, "input.json", fsize) catch break :default;

        try parse(filebuf); allocator.free(filebuf); 
    }}

    try std.io.getStdOut().writer().print("\nTOTAL EXECUTION TIME: {}\n", .{std.fmt.fmtDuration(timer.read())});

    // const A = try Matrix(f64).init(2, 2, allocator); defer A.deinit();
    //
    // A.randn(0, 1, 2);
    //
    // A.ptr(0, 0).* = 7; A.ptr(0, 1).* = 1;
    // A.ptr(1, 0).* = -1; A.ptr(1, 1).* = 3;
    //
    // var B = try Matrix(f64).init(A.rows, A.cols, allocator); defer B.deinit();
    // var C = try Matrix(f64).init(A.rows, A.cols, allocator); defer C.deinit();
    // var AJI = try Vector(f64).init(A.rows, allocator); defer AJI.deinit();
    // var AJR = try Vector(f64).init(A.rows, allocator); defer AJR.deinit();
    // var Q   = try Matrix(f64).init(A.rows, A.cols, allocator); defer Q.deinit();
    // var D   = try Matrix(f64).init(A.rows, A.cols, allocator); defer D.deinit();
    //
    // lapack.dgees(&Q, &D, A, &AJR, &AJI);
    //
    // blas.dgemm(&C, D, false, Q, true); blas.dgemm(&B, Q, false, C, false);
    //
    // try A.print(std.io.getStdOut().writer());
    // try Q.print(std.io.getStdOut().writer());
    // try D.print(std.io.getStdOut().writer());
    // try B.print(std.io.getStdOut().writer());
}
