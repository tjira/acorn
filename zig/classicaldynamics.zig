const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

pub fn ClassicalDynamics(comptime T: type) type {
    return struct {
        iterations: u32,
        mass: T,
        time_step: T,

        pub fn run(self: ClassicalDynamics(T), potential: fn (r: Vector(T)) anyerror!Matrix(T), r0: Matrix(T), v0: Matrix(T), a0: Matrix(T), s0: Vector(u8)) !void {
            if (r0.rows != v0.rows or r0.rows != a0.rows or r0.rows != s0.rows) return error.IncompatibleInitialConditions;
            for (0..r0.rows) |i| {
                try self.runTrajectory(potential, try r0.row(i), try v0.row(i), try a0.row(i), try s0.at(i));
            }
        }
        fn runTrajectory(self: ClassicalDynamics(T), potential: fn (r: Vector(T)) anyerror!Matrix(T), r0: Vector(T), v0: Vector(T), a0: Vector(T), s0: u8) !void {
            var r = try r0.clone(); var v = try v0.clone(); var a = try a0.clone(); const s = s0;
            for (0..self.iterations) |_| {
                const rp = try r.clone(); const vp = try v.clone(); const ap = try a.clone();// const sp = s;
                const dUs = try potentialDerivative(T, potential, r, s);
                a = try (try dUs.mulScalar(-1)).divScalar(self.mass);
                v = try vp.add(try (try ap.add(a)).mulScalar(0.5 * self.time_step));
                r = try rp.add(try (try v.add(try a.mulScalar(0.5 * self.time_step))).mulScalar(self.time_step));
                std.debug.print("{d:24.16}\n", .{try r.at(0)});
            }
        }
    };
}

fn potentialDerivative(comptime T: type, potential: fn (r: Vector(T)) anyerror!Matrix(T), r: Vector(T), s: u8) !Vector(T) {
    var dUs = try Vector(T).zero(r.rows, r.data.allocator);
    const spatial_step = 0.001;
    for (0..r.rows) |i| {
        const potential_plus = try potential(try r.addScalar(spatial_step));
        const potential_minus = try potential(try r.subScalar(spatial_step));
        try dUs.set(i, 0.5 * (try potential_plus.at(s, s) - try potential_minus.at(s, s)) / spatial_step);
    }
    return dUs;
}
