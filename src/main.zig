const std = @import("std"); const builtin = @import("builtin"); const Complex = std.math.complex.Complex;

const allocator = std.heap.page_allocator; const fsize = 2048;

pub const cdn = @import("classicaldynamics.zig"       );
pub const cim = @import("configurationinteraction.zig");
pub const ftr = @import("fouriertransform.zig"        );
pub const hfm = @import("hartreefock.zig"             );
pub const mat = @import("matrix.zig"                  );
pub const mpm = @import("mollerplesset.zig"           );
pub const mpt = @import("modelpotential.zig"          );
pub const qdn = @import("quantumdynamics.zig"         );
pub const vec = @import("vector.zig"                  );

pub const Matrix = @import("matrix.zig").Matrix;
pub const Tensor = @import("tensor.zig").Tensor;
pub const Vector = @import("vector.zig").Vector;

pub fn parse(filebuf: []const u8) !void {
    const inputjs = try std.json.parseFromSlice(std.json.Value, allocator, filebuf, .{}); defer inputjs.deinit();

    if (inputjs.value.object.contains("hartree_fock")) {

        const obj = try std.json.parseFromValue(hfm.HartreeFockOptions(f64), allocator, inputjs.value.object.get("hartree_fock").?, .{}); defer obj.deinit();

        const output = try hfm.run(f64, obj.value, true, allocator); defer output.deinit();
    }

    if (inputjs.value.object.contains("moller_plesset")) {

        const obj = try std.json.parseFromValue(mpm.MollerPlessetOptions(f64), allocator, inputjs.value.object.get("moller_plesset").?, .{}); defer obj.deinit();

        const output = try mpm.run(f64, obj.value, true, allocator); defer output.deinit();
    }

    if (inputjs.value.object.contains("configuration_interaction")) {

        const obj = try std.json.parseFromValue(cim.ConfigurationInteractionOptions(f64), allocator, inputjs.value.object.get("configuration_interaction").?, .{}); defer obj.deinit();

        const output = try cim.run(f64, obj.value, true, allocator); defer output.deinit();
    }

    if (inputjs.value.object.contains("classical_dynamics")) {

        const obj = try std.json.parseFromValue(cdn.ClassicalDynamicsOptions(f64), allocator, inputjs.value.object.get("classical_dynamics").?, .{}); defer obj.deinit();

        const output = try cdn.run(f64, obj.value, true, allocator); defer output.deinit();
    }

    if (inputjs.value.object.contains("quantum_dynamics")) {

        const obj = try std.json.parseFromValue(qdn.QuantumDynamicsOptions(f64), allocator, inputjs.value.object.get("quantum_dynamics").?, .{}); defer obj.deinit();

        const output = try qdn.run(f64, obj.value, true, allocator); defer output.deinit();
    }

    if (inputjs.value.object.contains("model_potential")) {

        const obj = try std.json.parseFromValue(mpt.ModelPotentialOptions(f64), allocator, inputjs.value.object.get("model_potential").?, .{}); defer obj.deinit();

        try mpt.write(f64, obj.value, allocator);
    }
}

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

        try std.io.getStdOut().writer().print("\nPROCESSED FILE: {s}\n", .{"input.json"});

        try parse(filebuf); allocator.free(filebuf); 
    }}

    std.debug.print("\nTOTAL EXECUTION TIME: {}\n", .{std.fmt.fmtDuration(timer.read())});
}
