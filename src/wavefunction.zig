const std = @import("std"); const Complex = std.math.Complex;

const ftr = @import("fouriertransform.zig");
const mat = @import("matrix.zig"          );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

pub fn Wavefunction(comptime T: type) type {
    return struct {
        data: Matrix(Complex(T)), shape: []usize, ndim: u32, nstate: u32, points: u32, allocator: std.mem.Allocator,

        pub fn init(ndim: u32, nstate: u32, points: u32, allocator: std.mem.Allocator) !Wavefunction(T) {
            const W = Wavefunction(T){
                .data = try Matrix(Complex(T)).init(std.math.pow(u32, points, ndim), nstate, allocator),
                .shape = try allocator.alloc(usize, ndim),
                .ndim = ndim,
                .nstate = nstate,
                .points = points,
                .allocator = allocator
            };

            for (0..ndim) |i| {
                W.shape[i] = points;
            }

            return W;
        }
        pub fn deinit(self: Wavefunction(T)) void {
            self.data.deinit(); self.allocator.free(self.shape);
        }
    };
}

pub fn adiabatize(comptime T: type, WA: *Wavefunction(T), W: Wavefunction(T), VC: std.ArrayList(Matrix(Complex(T)))) void {
    for (0..W.data.rows) |i| {

        var rowi =  W.data.row(i); rowi.rows = W.nstate; rowi.cols = 1;
        var rowo = WA.data.row(i); rowo.rows = W.nstate; rowo.cols = 1;

        mat.mam(Complex(T), &rowo, VC.items[i], rowi);
    }
}

pub fn density(comptime T: type, P: *Matrix(T), W: Wavefunction(T), dr: T) void {
    P.fill(0);

    for (0..W.nstate) |i| for (0..W.nstate) |j| for (0..W.data.rows) |k| {
        P.ptr(i, j).* += W.data.at(k, i).mul(W.data.at(k, j).conjugate()).re * dr;
    };
}

pub fn ekin(comptime T: type, W: Wavefunction(T), kvec: Matrix(T), mass: T, dr: T, T1: *Matrix(Complex(T)), T2: *Matrix(Complex(T))) T {
    var Ekin: T = 0;

    _ = T2;

    for (0..W.nstate) |i| {

        for (0..W.data.rows) |j| T1.ptr(j, 0).* = W.data.at(j, i);

        ftr.fftn(T, T1.data, W.shape, -1);

        for (0..W.data.rows) |j| {
            
            var ksqsum: T = 0;

            for (0..W.ndim) |k| ksqsum += kvec.at(j, k) * kvec.at(j, k);

            T1.ptr(j, 0).* = T1.at(j, 0).mul(Complex(T).init(ksqsum, 0));
        }

        ftr.fftn(T, T1.data, W.shape, 1);

        for (0..W.data.rows) |j| Ekin += T1.at(j, 0).mul(W.data.at(j, i).conjugate()).re * dr;
    }

    return 0.5 * Ekin / mass;
}

pub fn epot(comptime T: type, W: Wavefunction(T), V: std.ArrayList(Matrix(Complex(T))), dr: T) T {
    var Epot: T = 0;

    for (0..W.data.cols) |i| for (0..W.data.cols) |j| for (0..W.data.rows) |k| {
        Epot += W.data.at(k, i).conjugate().mul(V.items[k].at(i, j).mul(W.data.at(k, j))).re;
    };

    return Epot * dr;
}

pub fn guess(comptime T: type, W: *Wavefunction(T), rvec: Matrix(T), r: []const T, p: []const T, gamma: T, state: u32) void {
    W.data.fill(Complex(T).init(0, 0));

    for (0..W.data.rows) |i| for (0..W.data.cols) |j| if (j == state) {

        var exp: T = 0;

        for (0..W.ndim) |k| exp -= 0.5 * gamma * (rvec.at(i, k) - r[k]) * (rvec.at(i, k) - r[k]);

        W.data.ptr(i, j).* = Complex(T).init(std.math.exp(exp), 0);
    };

    for (0..W.data.rows) |i| for (0..W.data.cols) |j| {
        W.data.ptr(i, j).* = W.data.at(i, j).mul(std.math.complex.exp(Complex(T).init(0, rvec.at(i, 0) * p[0])));
    };
}

pub fn momentum(comptime T: type, p: *Vector(T), W: Wavefunction(T), kvec: Matrix(T), dr: T, T1: *Matrix(Complex(T)), T2: *Matrix(Complex(T))) void {
    p.fill(0);

    _ = T2;

    for (0..W.nstate) |i| {

        for (0..W.data.rows) |j| T1.ptr(j, 0).* = W.data.at(j, i);

        ftr.fftn(T, T1.data, W.shape, -1);

        for (0..W.data.rows) |j| T1.ptr(j, 0).* = T1.at(j, 0).mul(Complex(T).init(kvec.at(j, 0), 0));

        ftr.fftn(T, T1.data, W.shape, 1);

        for (0..W.data.rows) |j| p.ptr(0).* += T1.at(j, 0).mul(W.data.at(j, i).conjugate()).re * dr;
    }
}

pub fn normalize(comptime T: type, W: *Wavefunction(T), dr: T) void {
    const norm = std.math.sqrt(overlap(T, W.*, W.*, dr).magnitude());

    for (W.data.data) |*e| e.* = e.*.div(Complex(T).init(norm, 0));
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

    for (0..W.nstate) |i| for (0..W.data.rows) |j| for (0..rvec.cols) |k| {
        r.ptr(k).* += W.data.at(j, i).conjugate().mul(Complex(T).init(rvec.at(j, k), 0)).mul(W.data.at(j, i)).re * dr;
    };
}

pub fn propagate(comptime T: type, W: *Wavefunction(T), R: std.ArrayList(Matrix(Complex(T))), K: @TypeOf(R), T1: *Matrix(Complex(T)), T2: *Matrix(Complex(T))) void {
    _ = T2;

    for (0..W.data.rows) |i| {
        for (0..W.data.cols) |j| T1.data[j] = Complex(T).init(0, 0);
        for (0..W.data.cols) |j| for (0..W.data.cols) |k| {
            T1.data[j] = T1.data[j].add(R.items[i].at(j, k).mul(W.data.at(i, k)));
        };
        for (0..W.data.cols) |j| W.data.ptr(i, j).* = T1.data[j];
    }

    for (0..W.data.cols) |j| {
        for (0..W.data.rows) |i| T1.data[i] = W.data.at(i, j);
        ftr.fftn(T, T1.data, W.shape, -1);
        for (0..W.data.rows) |i| W.data.ptr(i, j).* = T1.at(i, 0);
    }

    for (0..W.data.rows) |i| {
        for (0..W.data.cols) |j| T1.data[j] = Complex(T).init(0, 0);
        for (0..W.data.cols) |j| for (0..W.data.cols) |k| {
            T1.data[j] = T1.data[j].add(K.items[i].at(j, k).mul(W.data.at(i, k)));
        };
        for (0..W.data.cols) |j| W.data.ptr(i, j).* = T1.data[j];
    }

    for (0..W.data.cols) |j| {
        for (0..W.data.rows) |i| T1.data[i] = W.data.at(i, j);
        ftr.fftn(T, T1.data, W.shape, 1);
        for (0..W.data.rows) |i| W.data.ptr(i, j).* = T1.at(i, 0);
    }

    for (0..W.data.rows) |i| {
        for (0..W.data.cols) |j| T1.data[j] = Complex(T).init(0, 0);
        for (0..W.data.cols) |j| for (0..W.data.cols) |k| {
            T1.data[j] = T1.data[j].add(R.items[i].at(j, k).mul(W.data.at(i, k)));
        };
        for (0..W.data.cols) |j| W.data.ptr(i, j).* = T1.data[j];
    }
}
