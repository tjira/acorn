const std = @import("std"); const Complex = std.math.Complex;

const ftr = @import("fouriertransform.zig");
const mat = @import("matrix.zig"          );

const Matrix = @import("matrix.zig").Matrix;

pub fn Wavefunction(comptime T: type) type {
    return struct {
        data: Matrix(Complex(T)), ndim: u32, nstate: u32, points: u32, allocator: std.mem.Allocator,

        pub fn init(ndim: u32, nstate: u32, points: u32, allocator: std.mem.Allocator) !Wavefunction(T) {
            return Wavefunction(T){.data = try Matrix(Complex(T)).init(std.math.pow(u32, points, ndim), nstate, allocator), .ndim = ndim, .nstate = nstate, .points = points, .allocator = allocator};
        }
        pub fn deinit(self: Wavefunction(T)) void {
            self.data.deinit();
        }

        pub fn normalize(self: Wavefunction(T), dr: T) void {
            const norm = std.math.sqrt(overlap(T, self, self, dr)); for (self.data.data) |*e| e.* = e.*.div(Complex(T).init(norm, 0));
        }
    };
}

pub fn density(comptime T: type, P: *Matrix(T), W: Wavefunction(T), dr: T) void {
    P.fill(0);

    for (0..W.nstate) |i| for (0..W.nstate) |j| for (0..W.data.rows) |k| {
        P.ptr(i, j).* += W.data.at(k, i).mul(W.data.at(k, j).conjugate()).re * dr;
    };
}

pub fn ekin(comptime T: type, W: Wavefunction(T), kvec: Matrix(T), mass: T, dr: T, T1: *Matrix(Complex(T)), T2: @TypeOf(T1)) T {
    var Ekin: T = 0;

    for (0..W.nstate) |i| {

        W.data.col(T1, i); ftr.dft(T, T2.data, T1.data, &[_]usize{W.data.rows}, -1);

        for (0..W.data.rows) |j| T2.ptr(j, 0).* = T2.at(j, 0).mul(Complex(T).init(kvec.at(j, 0) * kvec.at(j, 0), 0));

        ftr.dft(T, T1.data, T2.data, &[_]usize{W.data.rows}, 1);

        for (0..W.data.rows) |j| Ekin += T1.at(j, 0).mul(W.data.at(j, i).conjugate()).re;
    }

    return 0.5 * Ekin / mass * dr;
}

pub fn epot(comptime T: type, W: Wavefunction(T), V: std.ArrayList(Matrix(Complex(T))), dr: T) T {
    var Epot: T = 0;

    for (0..W.data.cols) |i| for (0..W.data.cols) |j| for (0..W.data.rows) |k| {
        Epot += W.data.at(k, j).conjugate().mul(V.items[k].at(i, j).mul(W.data.at(k, j))).re;
    };

    return Epot * dr;
}

pub fn guess(comptime T: type, W: *Wavefunction(T), rvec: Matrix(T), position: []const T, momentum: []const T, state: u32) void {
    W.data.fill(Complex(T).init(0, 0));

    for (0..W.data.rows) |i| for (0..W.data.cols) |j| if (j == state) {
        W.data.ptr(i, j).* = Complex(T).init(std.math.exp(-(rvec.at(i, 0) - position[0]) * (rvec.at(i, 0) - position[0])), 0);
    };

    for (0..W.data.rows) |i| for (0..W.data.cols) |j| {
        W.data.ptr(i, j).* = W.data.at(i, j).mul(std.math.complex.exp(Complex(T).init(0, rvec.at(i, 0) * momentum[0])));
    };
}

pub fn overlap(comptime T: type, W1: Wavefunction(T), W2: Wavefunction(T), dr: T) T {
    var s: T = 0;

    for (0..W1.nstate) |i| for (0..W1.nstate) |j| for (0..W1.data.rows) |k| {
        s += W1.data.at(k, i).conjugate().mul(W2.data.at(k, j)).re;
    };

    return s * dr;
}

pub fn propagate(comptime T: type, W: *Wavefunction(T), R: std.ArrayList(Matrix(Complex(T))), K: @TypeOf(R), T1: *Matrix(Complex(T)), T2: @TypeOf(T1)) void {
    var value = Complex(T).init(0, 0);

    for (0..W.data.rows) |i| for (0..W.data.cols) |j| {
        value.re = 0; value.im = 0; for (0..W.data.cols) |k| {value = value.add(R.items[i].at(j, k).mul(W.data.at(i, k)));} W.data.ptr(i, j).* = value;
    };

    for (0..W.data.cols) |j| {W.data.col(T1, j); ftr.dft(T, T2.data, T1.data, &[_]usize{W.data.rows}, 1); W.data.setcol(j, T2.*);}

    for (0..W.data.rows) |i| for (0..W.data.cols) |j| {
        value.re = 0; value.im = 0; for (0..W.data.cols) |k| {value = value.add(K.items[i].at(j, k).mul(W.data.at(i, k)));} W.data.ptr(i, j).* = value;
    };

    for (0..W.data.cols) |j| {W.data.col(T1, j); ftr.dft(T, T2.data, T1.data, &[_]usize{W.data.rows}, -1); W.data.setcol(j, T2.*);}

    for (0..W.data.rows) |i| for (0..W.data.cols) |j| {
        value.re = 0; value.im = 0; for (0..W.data.cols) |k| {value = value.add(R.items[i].at(j, k).mul(W.data.at(i, k)));} W.data.ptr(i, j).* = value;
    };
}
