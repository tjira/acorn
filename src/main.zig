const std = @import("std"); const builtin = @import("builtin");

const allocator = std.heap.page_allocator; const fsize = 2048;

const cdn = @import("classicaldynamics.zig");
const ftr = @import("fouriertransform.zig" );
const hfm = @import("hartreefock.zig"      );
const mat = @import("matrix.zig"           );
const mpt = @import("modelpotential.zig"   );
const qdn = @import("quantumdynamics.zig"  );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const CDO = @import("classicaldynamics.zig").ClassicalDynamicsOptions;
const HFO = @import("hartreefock.zig"      ).HartreeFockOptions      ;
const MPO = @import("modelpotential.zig"   ).ModelPotentialOptions   ;
const QDO = @import("quantumdynamics.zig"  ).QuantumDynamicsOptions  ;

pub fn main() !void {
    var timer = try std.time.Timer.start(); var args = try std.process.argsWithAllocator(allocator); defer args.deinit();

    try std.io.getStdOut().writer().print("ZIG VERSION: {}\n", .{builtin.zig_version});

    _ = args.next(); while (args.next()) |arg| {

        const filebuf = try std.fs.cwd().readFileAlloc(allocator, arg, fsize); defer allocator.free(filebuf);

        const inputjs = try std.json.parseFromSlice(std.json.Value, allocator, filebuf, .{}); defer inputjs.deinit();

        try std.io.getStdOut().writer().print("\nPROCESSED FILE: {s}\n", .{arg});

        if (inputjs.value.object.contains("restricted_hartree_fock")) {

            const obj = try std.json.parseFromValue(HFO(f64), allocator, inputjs.value.object.get("restricted_hartree_fock").?, .{}); defer obj.deinit();

            const output = try hfm.run(f64, obj.value, true, allocator); defer output.deinit();
        }

        if (inputjs.value.object.contains("classical_dynamics")) {

            const obj = try std.json.parseFromValue(CDO(f64), allocator, inputjs.value.object.get("classical_dynamics").?, .{}); defer obj.deinit();

            const output = try cdn.run(f64, obj.value, true, allocator); defer output.deinit();
        }

        if (inputjs.value.object.contains("quantum_dynamics")) {

            const obj = try std.json.parseFromValue(QDO(f64), allocator, inputjs.value.object.get("quantum_dynamics").?, .{}); defer obj.deinit();

            const output = try qdn.run(f64, obj.value, true, allocator); defer output.deinit();
        }

        if (inputjs.value.object.contains("model_potential")) {

            const obj = try std.json.parseFromValue(MPO(f64), allocator, inputjs.value.object.get("model_potential").?, .{}); defer obj.deinit();

            try mpt.write(f64, obj.value, allocator);
        }
    }

    std.debug.print("\nTOTAL EXECUTION TIME: {}\n", .{std.fmt.fmtDuration(timer.read())});
}

test {
    std.testing.refAllDecls(@This());
}
