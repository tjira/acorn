//! Module for the wavefunction object.

const std = @import("std"); const cwp = @import("cwrapper.zig");

const ftr = @import("fouriertransform.zig");
const mat = @import("matrix.zig"          );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

/// The wavefunction object.
pub fn Wavefunction(comptime T: type) type {
    return struct {
        data: Matrix(std.math.Complex(T)), shape: []usize, ndim: usize, npoint: usize, nstate: usize, dr: T, allocator: std.mem.Allocator,

        /// Initialize the wavefunction object with the given number of dimensions, states, and points.
        pub fn init(ndim: usize, nstate: usize, points: usize, dr: T, allocator: std.mem.Allocator) !Wavefunction(T) {
            const W = Wavefunction(T){
                .data = try Matrix(std.math.Complex(T)).init(std.math.pow(usize, points, ndim), nstate, allocator),
                .shape = try allocator.alloc(usize, ndim),
                .ndim = ndim,
                .npoint = std.math.pow(usize, points, ndim),
                .nstate = nstate,
                .dr = dr,
                .allocator = allocator
            };

            for (0..ndim) |i| {
                W.shape[i] = @intCast(points);
            }

            return W;
        }

        /// Free the memory allocated for the wavefunction object.
        pub fn deinit(self: Wavefunction(T)) void {
            self.data.deinit(); self.allocator.free(self.shape);
        }

        /// Index the underlying data matrix.
        pub fn at(self: Wavefunction(T), i: usize, j: usize) std.math.Complex(T) {
            return self.data.at(i, j);
        }

        /// Copy the self data to the destination wavefunction.
        pub fn memcpy(self: Wavefunction(T), dest: Wavefunction(T)) void {
            self.data.memcpy(dest.data);
        }

        /// Normalize the wavefunction.
        pub fn normalize(self: *Wavefunction(T)) void {
            mat.divs(std.math.Complex(T), &self.data, self.data, std.math.Complex(T).init(std.math.sqrt(self.overlap(self.*).magnitude()), 0));
        }

        /// Calculate the overlap with another wavefunction.
        pub fn overlap(self: Wavefunction(T), W: Wavefunction(T)) std.math.Complex(T) {
            var s = std.math.Complex(T).init(0, 0);

            for (0..self.nstate) |i| for (0..self.nstate) |j| for (0..self.npoint) |k| {
                s = s.add(self.at(k, i).conjugate().mul(W.at(k, j)));
            };

            return s.mul(std.math.Complex(T).init(self.dr, 0));
        }

        /// Get a pointer to the data at the given row and column.
        pub fn ptr(self: Wavefunction(T), i: usize, j: usize) *std.math.Complex(T) {
            return self.data.ptr(i, j);
        }

        /// Returns a row of the wavefunction as a matrix.
        pub fn row(self: Wavefunction(T), i: usize) Matrix(std.math.Complex(T)) {
            return self.data.row(i);
        }
    };
}

/// Given the transformation matrices for each point in the grid VC, this function adiabatizes the wavefunction W and stores the result in WA.
pub fn adiabatize(comptime T: type, WA: *Wavefunction(T), W: Wavefunction(T), VC: std.ArrayList(Matrix(std.math.Complex(T)))) !void {
    for (0..W.npoint) |i| {
        var rowa = WA.row(i).vector().matrix(); try cwp.Blas(T).zgemm(&rowa, VC.items[i], true, W.row(i).vector().matrix(), false);
    }
}

/// Calculate the density matrix of the wavefunction W and store the result in P.
pub fn density(comptime T: type, P: *Matrix(std.math.Complex(T)), W: Wavefunction(T)) void {
    P.fill(std.math.Complex(T).init(0, 0));

    for (0..W.nstate) |i| for (0..W.nstate) |j| for (0..W.npoint) |k| {
        P.ptr(i, j).* = P.at(i, j).add(W.at(k, i).mul(W.at(k, j).conjugate()).mul(std.math.Complex(T).init(W.dr, 0)));
    };
}

/// Given the transformation matrices for each point in the grid VC, this function adiabatizes the wavefunction W and stores the result in WA.
pub fn diabatize(comptime T: type, WD: *Wavefunction(T), W: Wavefunction(T), VC: std.ArrayList(Matrix(std.math.Complex(T)))) !void {
    for (0..W.npoint) |i| {
        var rowd = WD.row(i).vector().matrix(); try cwp.Blas(T).zgemm(&rowd, VC.items[i], false, W.row(i).vector().matrix(), false);
    }
}

/// Calculate the kinetic energy of the wavefunction W. This function needs the grid in the k-space kvec.
pub fn ekin(comptime T: type, W: Wavefunction(T), kvec: Matrix(T), mass: T, T1: *Matrix(std.math.Complex(T))) !T {
    var Ekin: T = 0;

    for (0..W.nstate) |i| {

        for (0..W.npoint) |j| T1.ptr(j, 0).* = W.at(j, i);

        try ftr.fftn(f64, T1.data, W.shape, -1);

        for (0..W.npoint) |j| {
            
            var ksqsum: T = 0;

            for (0..W.ndim) |k| ksqsum += kvec.at(j, k) * kvec.at(j, k);

            T1.ptr(j, 0).* = T1.at(j, 0).mul(std.math.Complex(T).init(ksqsum, 0));
        }

        try ftr.fftn(f64, T1.data, W.shape, 1);

        for (0..W.npoint) |j| Ekin += T1.at(j, 0).mul(W.at(j, i).conjugate()).re * W.dr;
    }

    return 0.5 * Ekin / mass;
}

/// Calculate the potential energy of the wavefunction W. This function needs the potential matrices at each point in the grid V.
pub fn epot(comptime T: type, W: Wavefunction(T), V: std.ArrayList(Matrix(std.math.Complex(T)))) T {
    var Epot: T = 0;

    for (0..W.nstate) |i| for (0..W.nstate) |j| for (0..W.npoint) |k| {
        Epot += W.at(k, i).conjugate().mul(V.items[k].at(i, j).mul(W.at(k, j))).re;
    };

    return Epot * W.dr;
}

/// Calculate the guess wavefunction as a gaussian centered at r with momentum p on a given state.
pub fn guess(comptime T: type, W: *Wavefunction(T), rvec: Matrix(T), r: []const T, p: []const T, gamma: T, state: usize) void {
    W.data.fill(std.math.Complex(T).init(0, 0));

    for (0..W.npoint) |i| for (0..W.nstate) |j| if (j == state) {

        var exp: T = 0;

        for (0..W.ndim) |k| exp -= 0.5 * gamma * (rvec.at(i, k) - r[k]) * (rvec.at(i, k) - r[k]);

        W.ptr(i, j).* = std.math.Complex(T).init(std.math.exp(exp), 0);
    };

    for (0..W.npoint) |i| for (0..W.nstate) |j| {

        var exp: T = 0;

        for (0..W.ndim) |k| exp += rvec.at(i, k) * p[k];

        W.ptr(i, j).* = W.at(i, j).mul(std.math.complex.exp(std.math.Complex(T).init(0, exp)));
    };
}

/// Calculate the momentum of the wavefunction W. This function needs the grid in the k-space kvec.
pub fn momentum(comptime T: type, p: *Vector(T), W: Wavefunction(T), kvec: Matrix(T), T1: *Matrix(std.math.Complex(T))) !void {
    p.fill(0);

    for (0..W.nstate) |i| {

        for (0..kvec.cols) |k| {

            for (0..W.npoint) |j| T1.ptr(j, 0).* = W.at(j, i);

            try ftr.fftn(f64, T1.data, W.shape, -1);

            for (0..W.npoint) |j| T1.ptr(j, 0).* = T1.at(j, 0).mul(std.math.Complex(T).init(kvec.at(j, k), 0));

            try ftr.fftn(f64, T1.data, W.shape, 1);

            for (0..W.npoint) |j| p.ptr(k).* += T1.at(j, 0).mul(W.at(j, i).conjugate()).re * W.dr;
        }
    }
}

/// Calculate the position of the wavefunction W. This function needs the grid in the r-space rvec.
pub fn position(comptime T: type, r: *Vector(T), W: Wavefunction(T), rvec: Matrix(T)) void {
    r.fill(0);

    for (0..W.nstate) |i| for (0..W.npoint) |j| for (0..rvec.cols) |k| {
        r.ptr(k).* += W.at(j, i).conjugate().mul(std.math.Complex(T).init(rvec.at(j, k), 0)).mul(W.at(j, i)).re * W.dr;
    };
}

/// Propagate the wavefunction W using the propagation matrices for each point of the grid R and K.
pub fn propagate(comptime T: type, W: *Wavefunction(T), R: std.ArrayList(Matrix(std.math.Complex(T))), K: @TypeOf(R), T1: *Matrix(std.math.Complex(T))) !void {
    for (0..W.npoint) |i| {
        for (0..W.nstate) |j| T1.ptr(j, 0).* = std.math.Complex(T).init(0, 0);
        for (0..W.nstate) |j| for (0..W.nstate) |k| {
            T1.ptr(j, 0).* = T1.at(j, 0).add(R.items[i].at(j, k).mul(W.at(i, k)));
        };
        for (0..W.nstate) |j| W.ptr(i, j).* = T1.at(j, 0);
    }

    for (0..W.nstate) |j| {
        for (0..W.npoint) |i| T1.data[i] = W.at(i, j);
        try ftr.fftn(f64, T1.data, W.shape, -1);
        for (0..W.npoint) |i| W.ptr(i, j).* = T1.at(i, 0);
    }

    for (0..W.npoint) |i| {
        for (0..W.nstate) |j| T1.ptr(j, 0).* = std.math.Complex(T).init(0, 0);
        for (0..W.nstate) |j| for (0..W.nstate) |k| {
            T1.ptr(j, 0).* = T1.at(j, 0).add(K.items[i].at(j, k).mul(W.at(i, k)));
        };
        for (0..W.nstate) |j| W.ptr(i, j).* = T1.at(j, 0);
    }

    for (0..W.nstate) |j| {
        for (0..W.npoint) |i| T1.data[i] = W.at(i, j);
        try ftr.fftn(f64, T1.data, W.shape, 1);
        for (0..W.npoint) |i| W.ptr(i, j).* = T1.at(i, 0);
    }

    for (0..W.npoint) |i| {
        for (0..W.nstate) |j| T1.ptr(j, 0).* = std.math.Complex(T).init(0, 0);
        for (0..W.nstate) |j| for (0..W.nstate) |k| {
            T1.ptr(j, 0).* = T1.at(j, 0).add(R.items[i].at(j, k).mul(W.at(i, k)));
        };
        for (0..W.nstate) |j| W.ptr(i, j).* = T1.at(j, 0);
    }
}
