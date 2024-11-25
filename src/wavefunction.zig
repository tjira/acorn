const std = @import("std"); const Complex = std.math.Complex; const gsl_fft = @cImport(@cInclude("gsl/gsl_fft_complex.h"));

const ftr = @import("fouriertransform.zig");
const mat = @import("matrix.zig"          );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

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
            const norm = std.math.sqrt(overlap(T, self, self, dr).magnitude()); for (self.data.data) |*e| e.* = e.*.div(Complex(T).init(norm, 0));
        }
    };
}

pub fn adiabatize(comptime T: type, WA: *Wavefunction(T), W: Wavefunction(T), VC: std.ArrayList(Matrix(Complex(T)))) void {
    WA.data.fill(Complex(T).init(0, 0));

    for (0..W.data.rows) |i| for (0..W.nstate) |j| for (0..W.nstate) |k| {
        WA.data.ptr(i, j).* = WA.data.at(i, j).add(VC.items[i].at(k, j).conjugate().mul(W.data.at(i, k)));
    };
}

pub fn density(comptime T: type, P: *Matrix(T), W: Wavefunction(T), dr: T) void {
    for (0..W.nstate) |i| for (0..W.nstate) |j| {
        var pij = Complex(T).init(0, 0); for (0..W.data.rows) |k| {pij = pij.add(W.data.at(k, i).mul(W.data.at(k, j).conjugate()));} P.ptr(i, j).* = pij.magnitude() * dr;
    };
}

pub fn ekin(comptime T: type, W: Wavefunction(T), kvec: Matrix(T), mass: T, dr: T, WS: *Matrix(Complex(T)), GSLFT: *gsl_fft.gsl_fft_complex_wavetable, GSLFW: *gsl_fft.gsl_fft_complex_workspace) T {
    var Ekin: T = 0;

    for (0..W.nstate) |i| {

        for (0..W.data.rows) |j| {WS.ptr(j, 0).* = W.data.at(j, i);} ftr.fft(T, WS.data, -1, GSLFT, GSLFW);

        for (0..W.data.rows) |j| WS.ptr(j, 0).* = WS.at(j, 0).mul(Complex(T).init(kvec.at(j, 0) * kvec.at(j, 0), 0));

        ftr.fft(T, WS.data, 1, GSLFT, GSLFW);

        for (0..W.data.rows) |j| Ekin += WS.at(j, 0).mul(W.data.at(j, i).conjugate()).re;
    }

    return 0.5 * Ekin / mass * dr;
}

pub fn epot(comptime T: type, W: Wavefunction(T), V: std.ArrayList(Matrix(Complex(T))), dr: T) T {
    var Epot: T = 0;

    for (0..W.data.cols) |i| for (0..W.data.cols) |j| for (0..W.data.rows) |k| {
        Epot += W.data.at(k, i).conjugate().mul(V.items[k].at(i, j).mul(W.data.at(k, j))).re;
    };

    return Epot * dr;
}

pub fn guess(comptime T: type, W: *Wavefunction(T), rvec: Matrix(T), r: []const T, p: []const T, state: u32) void {
    W.data.fill(Complex(T).init(0, 0));

    for (0..W.data.rows) |i| for (0..W.data.cols) |j| if (j == state) {
        W.data.ptr(i, j).* = Complex(T).init(std.math.exp(-(rvec.at(i, 0) - r[0]) * (rvec.at(i, 0) - r[0])), 0);
    };

    for (0..W.data.rows) |i| for (0..W.data.cols) |j| {
        W.data.ptr(i, j).* = W.data.at(i, j).mul(std.math.complex.exp(Complex(T).init(0, rvec.at(i, 0) * p[0])));
    };
}

pub fn momentum(comptime T: type, p: *Vector(T), W: Wavefunction(T), kvec: Matrix(T), dr: T, WS: *Matrix(Complex(T)), GSLFT: *gsl_fft.gsl_fft_complex_wavetable, GSLFW: *gsl_fft.gsl_fft_complex_workspace) void {
    p.fill(0);

    for (0..W.nstate) |i| {

        for (0..W.data.rows) |j| {WS.ptr(j, 0).* = W.data.at(j, i);} ftr.fft(T, WS.data, -1, GSLFT, GSLFW);

        for (0..W.data.rows) |j| WS.ptr(j, 0).* = WS.at(j, 0).mul(Complex(T).init(kvec.at(j, 0), 0));

        ftr.fft(T, WS.data, 1, GSLFT, GSLFW);

        for (0..W.data.rows) |j| p.ptr(0).* += WS.at(j, 0).mul(W.data.at(j, i).conjugate()).re * dr;
    }
}

pub fn overlap(comptime T: type, W1: Wavefunction(T), W2: Wavefunction(T), dr: T) Complex(T) {
    var s = Complex(T).init(0, 0);

    for (0..W1.nstate) |i| for (0..W1.nstate) |j| for (0..W1.data.rows) |k| {
        s = s.add(W1.data.at(k, i).conjugate().mul(W2.data.at(k, j)));
    };

    return s.mul(Complex(T).init(dr, 0));
}

pub fn position(comptime T: type, r: *Vector(T), W: Wavefunction(T), rvec: Matrix(T), dr: T) void {
    r.fill(0);

    for (0..W.nstate) |i| for (0..W.data.rows) |j| {
        r.ptr(0).* += W.data.at(j, i).conjugate().mul(Complex(T).init(rvec.at(j, 0), 0)).mul(W.data.at(j, i)).re * dr;
    };
}

pub fn propagate(comptime T: type, W: *Wavefunction(T), R: std.ArrayList(Matrix(Complex(T))), K: @TypeOf(R), WS: *Matrix(Complex(T)), GSLFT: *gsl_fft.gsl_fft_complex_wavetable, GSLFW: *gsl_fft.gsl_fft_complex_workspace) void {
    for (0..W.data.rows) |i| {
        for (0..W.data.cols) |j| WS.data[j] = Complex(T).init(0, 0);
        for (0..W.data.cols) |j| for (0..W.data.cols) |k| {WS.data[j] = WS.data[j].add(R.items[i].at(j, k).mul(W.data.at(i, k)));};
        for (0..W.data.cols) |j| W.data.ptr(i, j).* = WS.data[j];
    }

    for (0..W.data.cols) |j| {
        for (0..W.data.rows) |i| {WS.data[i] = W.data.at(i, j);} ftr.fft(T, WS.data, -1, GSLFT, GSLFW); for (0..W.data.rows) |i| {W.data.ptr(i, j).* = WS.at(i, 0);}
    }

    for (0..W.data.rows) |i| {
        for (0..W.data.cols) |j| WS.data[j] = Complex(T).init(0, 0);
        for (0..W.data.cols) |j| for (0..W.data.cols) |k| {WS.data[j] = WS.data[j].add(K.items[i].at(j, k).mul(W.data.at(i, k)));};
        for (0..W.data.cols) |j| W.data.ptr(i, j).* = WS.data[j];
    }

    for (0..W.data.cols) |j| {
        for (0..W.data.rows) |i| {WS.data[i] = W.data.at(i, j);} ftr.fft(T, WS.data, 1, GSLFT, GSLFW); for (0..W.data.rows) |i| {W.data.ptr(i, j).* = WS.at(i, 0);}
    }

    for (0..W.data.rows) |i| {
        for (0..W.data.cols) |j| WS.data[j] = Complex(T).init(0, 0);
        for (0..W.data.cols) |j| for (0..W.data.cols) |k| {WS.data[j] = WS.data[j].add(R.items[i].at(j, k).mul(W.data.at(i, k)));};
        for (0..W.data.cols) |j| W.data.ptr(i, j).* = WS.data[j];
    }
}
