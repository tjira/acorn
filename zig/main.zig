const std = @import("std");

const cdn = @import("classicaldynamics.zig");
const ftr = @import("fouriertransform.zig");
const mat = @import("matrix.zig"          );
const mpt = @import("modelpotential.zig"  );
const qdn = @import("quantumdynamics.zig"  );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const allocator = std.heap.page_allocator;

pub fn main() !void {
    var timer = try std.time.Timer.start();

    const path = "example/input/cdyn_tripleState1D-1.json";
    // const path = "example/input/qrdn_tripleState1D-1.json";
    // const path = "example/input/qidn_harmonic1D-1.json";

    const mode = cdn.ClassicalDynamicsOptions(f64);
    // const mode = qdn.QuantumDynamicsOptions(f64);

    const buffer = try std.fs.cwd().readFileAlloc(allocator, path, 2048); defer allocator.free(buffer);
    const parsed = try std.json.parseFromSlice(mode, allocator, buffer, .{}); defer parsed.deinit();
    const opt = parsed.value; std.debug.print("{s}\n", .{buffer});

    if (@TypeOf(opt) == cdn.ClassicalDynamicsOptions(f64)) try cdn.run(f64, opt, allocator);
    if (@TypeOf(opt) == qdn.QuantumDynamicsOptions(f64)  ) try qdn.run(f64, opt, allocator);

    const pot_start: f64 = -16; const pot_end: f64 = 16; const pot_points: f64 = 1024;
    try mpt.write(f64, "POTENTIAL-DIABATIC.mat",  opt.potential, pot_start, pot_end, pot_points, false, allocator);
    try mpt.write(f64, "POTENTIAL-ADIABATIC.mat", opt.potential, pot_start, pot_end, pot_points, true,  allocator);

    std.debug.print("\nTOTAL EXECUTION TIME: {}\n\n\n\n", .{std.fmt.fmtDuration(timer.read())});
}
