const std = @import("std");

const allocator = std.heap.page_allocator; const fsize = 2048;

const cdn = @import("classicaldynamics.zig");
const mat = @import("matrix.zig"           );
const mpt = @import("modelpotential.zig"   );
const qdn = @import("quantumdynamics.zig"  );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const CDO = @import("classicaldynamics.zig").ClassicalDynamicsOptions;
const QDO = @import("quantumdynamics.zig"  ).QuantumDynamicsOptions  ;
const MPO = @import("modelpotential.zig"   ).ModelPotentialOptions   ;

pub fn main() !void {
    var timer = try std.time.Timer.start(); var args = try std.process.argsWithAllocator(allocator); defer args.deinit();

    _ = args.next(); while (args.next()) |arg| {

        const filebuf = try std.fs.cwd().readFileAlloc(allocator, arg, fsize); defer allocator.free(filebuf);

        const inputjs = try std.json.parseFromSlice(std.json.Value, allocator, filebuf, .{}); defer inputjs.deinit();

        if (inputjs.value.object.contains("classical_dynamics")) {

            const obj = try std.json.parseFromValue(CDO(f64), allocator, inputjs.value.object.get("classical_dynamics").?, .{}); defer obj.deinit();

            try cdn.run(f64, obj.value, allocator);
        }

        if (inputjs.value.object.contains("quantum_dynamics")) {

            const obj = try std.json.parseFromValue(QDO(f64), allocator, inputjs.value.object.get("quantum_dynamics").?, .{}); defer obj.deinit();

            try qdn.run(f64, obj.value, allocator);
        }

        if (inputjs.value.object.contains("model_potential")) {

            const obj = try std.json.parseFromValue(MPO(f64), allocator, inputjs.value.object.get("model_potential").?, .{}); defer obj.deinit();

            try mpt.write(f64, obj.value, allocator);
        }
    }

    std.debug.print("\nTOTAL EXECUTION TIME: {}\n", .{std.fmt.fmtDuration(timer.read())});

    // const A = try Matrix(f64).init(2, 2, allocator); defer  A.deinit();
    // var  AC = try Matrix(f64).init(2, 2, allocator); defer AC.deinit();
    // var  AJ = try Matrix(f64).init(2, 2, allocator); defer AJ.deinit();
    //
    // var v = A.rowptr(0);
    // A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2; A.ptr(1, 0).* = 2; A.ptr(1, 1).* = 1;

    // v.ptr(0, 0).* = 12;

    // try A.print(std.io.getStdOut().writer());

    // mat.eigh(f64, &AJ, &AC, A);
    //
    // try A.print(std.io.getStdOut().writer());
    // std.debug.print("\n", .{});
    // try AJ.print(std.io.getStdOut().writer());
    // std.debug.print("\n", .{});
    // try AC.print(std.io.getStdOut().writer());
}
