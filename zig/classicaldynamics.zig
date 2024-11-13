const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const allocator = std.heap.page_allocator;

pub const ClassicalDynamics = struct {
    iterations: u32,
    mass: f64,
    timeStep: f64,

    pub fn run(self: ClassicalDynamics, potential: fn (r: Vector(f64)) anyerror!Matrix(f64), r0: Matrix(f64), v0: Matrix(f64), a0: Matrix(f64), s0: Vector(u8)) !void {
        if (r0.rows != v0.rows or r0.rows != a0.rows or r0.rows != s0.rows) return error.IncompatibleInitialConditions; // check for the same number of initial conditions
        for (0..r0.rows) |i| { // loop over all trajectories
            try self.runTrajectory(potential, try r0.row(i), try v0.row(i), try a0.row(i), try s0.at(i)); // run the trajectory
        }
    }
    fn runTrajectory(self: ClassicalDynamics, potential: fn (r: Vector(f64)) anyerror!Matrix(f64), r0: Vector(f64), v0: Vector(f64), a0: Vector(f64), s0: u8) !void {
        var r = try r0.clone(); var v = try v0.clone(); var a = try a0.clone(); const s = s0; // clone the initial position, velocity, acceleration and state
        for (0..self.iterations) |_| { // loop over the iteration count for propagation
            const rp = try r.clone(); const vp = try v.clone(); const ap = try a.clone();// const sp = s; // save the previous position, velocity, acceleration and state
            const dUs = try potentialDerivative(potential, r, s); // calculate the potential derivative at the current point and state
            a = try (try dUs.negate()).divScalar(self.mass); // calculate the acceleration at the current point
            v = try vp.add(try (try ap.add(a)).mulScalar(0.5 * self.timeStep)); // propagate the velocity
            r = try rp.add(try (try v.add(try a.mulScalar(0.5 * self.timeStep))).mulScalar(self.timeStep)); // propagate the position
            r.debug();
        }
    }
};

fn potentialDerivative(potential: fn (r: Vector(f64)) anyerror!Matrix(f64), r: Vector(f64), s: u8) !Vector(f64) {
    var dUs = try Vector(f64).zero(r.rows, allocator); // define the potential derivative for the current state
    const spatialStep = 0.001; // define the step used for the numerical derivative
    for (0..r.rows) |i| { // loop over all coordinates
        const potentialPlus = try potential(try r.addScalar(spatialStep)); // define the forward step
        const potentialMinus = try potential(try r.subScalar(spatialStep)); // define the forward step
        try dUs.set(i, 0.5 * (try potentialPlus.at(s, s) - try potentialMinus.at(s, s)) / spatialStep); // calculate and set the derivative
    }
    return dUs; // return the potential
}
