const std = @import("std");

const SurfaceHopping = @import("surfacehopping.zig").SurfaceHopping;
const Matrix         = @import("matrix.zig"        ).Matrix        ;
const Vector         = @import("vector.zig"        ).Vector        ;

pub fn ClassicalDynamics(comptime T: type) type {
    return struct {
        iterations: u32,
        mass: T,
        time_step: T,
        adiabatic: bool,
        seed: usize,

        // =============================================================================================================================================================================================

        pub const IterationData = struct {
            t: T,
            r: Vector(T),
            v: Vector(T),
            a: Vector(T),
            s: u8,
            Ekin: T,
            Epot: T,
            U: Matrix(T)
        };

        // =============================================================================================================================================================================================

        pub fn run(self: ClassicalDynamics(T), potential: fn (r: Vector(T)) anyerror!Matrix(T), r0: Matrix(T), v0: Matrix(T), a0: Matrix(T), s0: Vector(u8)) !void {
            if (r0.rows != v0.rows or r0.rows != a0.rows or r0.rows != s0.rows) return error.IncompatibleInitialConditions;
            for (0..r0.rows) |i| {
                const traj_data = try self.runTrajectory(potential, try r0.row(i), try v0.row(i), try a0.row(i), s0.at(i), i + 1);
                _ = traj_data;
            }
        }
        fn runTrajectory(self: ClassicalDynamics(T), potential: fn (r: Vector(T)) anyerror!Matrix(T), r0: Vector(T), v0: Vector(T), a0: Vector(T), s0: u8, index: usize) !std.ArrayList(IterationData) {
            var prng = std.Random.DefaultPrng.init(self.seed); const rand = prng.random();
            var r = try r0.clone(); var v = try v0.clone(); var a = try a0.clone(); var s = s0;
            var traj_data = std.ArrayList(IterationData).init(r.data.allocator);
            for (0..self.iterations) |i| {
                const rp = try r.clone(); const vp = try v.clone(); const ap = try a.clone(); const sp = s;
                const F = try (try potentialDerivative(T, potential, r, s)).mulScalar(-1);
                a = try F.divScalar(self.mass);
                v = try vp.add(try (try ap.add(a)).mulScalar(0.5 * self.time_step));
                r = try rp.add(try (try v.add(try a.mulScalar(0.5 * self.time_step))).mulScalar(self.time_step));
                const U = try potential(r);
                const Ekin = 0.5 * self.mass * (try v.mul(v)).sum();
                const Epot = U.at(s, s);
                const p = try SurfaceHopping(T).landauZener(traj_data, self.time_step, self.adiabatic);
                s = performStateJump(T, U, s, p, rand);
                v = rescaleVelocity(T, v, Ekin, U.at(s, s) - U.at(sp, sp));
                try traj_data.append(.{.t = @as(T, @floatFromInt(i + 1)) * self.time_step, .r = r, .v = v, .a = a, .s = s, .Ekin = Ekin, .Epot = Epot, .U = U});
                std.debug.print("{d:6} {d:6} {d:12.6} {d:12.6} {d:12.6} {d:1} {d:12.6}\n", .{index, i + 1, Ekin, Epot, Ekin + Epot, s, r.at(0)});
            }
            return traj_data;
        }
    };
}

fn performStateJump(comptime T: type, U: Matrix(T), s: u8, p: ?std.AutoHashMap(u8, T), rand: std.Random) u8 {
    _ = U;
    _ = rand;
    _ = p;
    return s;
}

fn potentialDerivative(comptime T: type, potential: fn (r: Vector(T)) anyerror!Matrix(T), r: Vector(T), s: u8) !Vector(T) {
    var dUs = try Vector(T).zero(r.rows, r.data.allocator);
    const spatial_step = 0.001;
    for (0..r.rows) |i| {
        const potential_plus = try potential(try r.addScalar(spatial_step));
        const potential_minus = try potential(try r.subScalar(spatial_step));
        dUs.set(i, 0.5 * (potential_plus.at(s, s) - potential_minus.at(s, s)) / spatial_step);
    }
    return dUs;
}

fn rescaleVelocity(comptime T: type, v: Vector(T), Ediff: T, Ekin: T) Vector(T) {
    _ = Ediff;
    _ = Ekin;
    return v;
}
