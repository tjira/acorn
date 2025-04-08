//! This module provides some specific linear algebra functions. It includes functions for finding eigenvalues and eigenvectors, solving linear systems, and performing QR decomposition.

const std = @import("std");

const mat = @import("matrix.zig");
const mth = @import("math.zig"  );
const vec = @import("vector.zig");

const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;



/// The Davidson algorithm for finding the eigenvalues and eigenvectors of a real symmetric matrix A.
pub fn davidson(comptime T: type, J: *Matrix(T), C: *Matrix(T), A: Matrix(T), k: usize, m: usize, allocator: std.mem.Allocator) !void {
    const tolerance: T = 1e-12;

    var JP = try Matrix(T).init(J.rows, J.cols, allocator); defer JP.deinit();
    var d  = try Vector(T).init(A.rows,         allocator); defer  d.deinit();

    var basis = std.ArrayList(Vector(T)).init(allocator); defer basis.deinit(); mat.diag(T, &d, A);

    for (0..k) |i| {
        var e = try Vector(T).init(A.rows, allocator); e.fill(0); e.ptr(i).* = 1; try basis.append(e);
    }

    for (0..J.cols) |i| J.ptr(i, i).* = std.math.inf(T);

    while (basis.items.len <= m) {
        var V  = try Matrix(T).init(A.rows,          basis.items.len, allocator); defer  V.deinit();
        var K  = try Matrix(T).init(basis.items.len, basis.items.len, allocator); defer  K.deinit();
        var KJ = try Matrix(T).init(basis.items.len, basis.items.len, allocator); defer KJ.deinit();
        var KC = try Matrix(T).init(basis.items.len, basis.items.len, allocator); defer KC.deinit();
        var T1 = try Matrix(T).init(basis.items.len, basis.items.len, allocator); defer T1.deinit();
        var T2 = try Matrix(T).init(basis.items.len, basis.items.len, allocator); defer T2.deinit();
        var T3 = try Matrix(T).init(basis.items.len, A.rows,          allocator); defer T3.deinit();

        for (0..basis.items.len) |i| for (0..A.rows) |j| {V.ptr(j, i).* = basis.items[i].at(j);};

        mat.mam(T, &T3, V, A); mat.mm(T, &K, T3, V); @memcpy(JP.data , J.data);

        try eighQr(T, &KJ, &KC, K, allocator);

        for (0..KJ.rows) |i| {J.ptr(i, i).* = KJ.at(i, i);} for (0..KC.rows) |i| for (0..KC.cols) |j| {C.ptr(i, j).* = KC.at(i, j);};

        var sumsq: T = 0; for (0..k) |i| sumsq += (J.at(i, i) - JP.at(i, i)) * (J.at(i, i) - JP.at(i, i));

        if (std.math.sqrt(sumsq) < tolerance) break;

        for (0..k) |i| {

            var x = try Vector(T).init(V.rows, allocator); x.fill(0);
            var r = try Vector(T).init(V.rows, allocator); r.fill(0);

            for (0..A.rows) |j| for (0..basis.items.len) |l| {
                x.ptr(j).* += V.at(j, l) * KC.at(l, i);
            };

            for (0..A.rows) |j| for (0..x.rows) |l| {
                r.ptr(j).* += A.at(j, l) * x.at(l);
            };

            for (0..x.rows) |j| {
                r.ptr(j).* -= J.at(i, i) * x.at(j); r.ptr(j).* /= (J.at(i, i) - d.at(j) + 1e-14);
            }

            for (0..basis.items.len) |j| {

                var dot: T = 0;

                for (0..A.rows) |l| {
                    dot += basis.items[j].at(l) * r.at(l);
                }

                for (0..A.rows) |l| {
                    r.ptr(l).* -= dot * basis.items[j].at(l);
                }
            }

            const rnorm = r.norm(); x.deinit();

            if (rnorm > 1e-14) {
                for (0..r.rows) |j| {r.ptr(j).* /= rnorm;} try basis.append(r);
            }
        }
    }

    for (basis.items) |v| v.deinit();
}

/// Find the eigenvalues and eigenvectors of a real symmetric matrix A. The eigenvalues are stored in the diagonal of the matrix J, and the eigenvectors are stored in the columns of the matrix C. The matrices T1 and T2 are temporary matrices used in the computation.
pub fn eighJacobi(comptime T: type, J: *Matrix(T), C: *Matrix(T), A: Matrix(T), T1: *Matrix(T), T2: *Matrix(T)) void {
    const tolerance: T = 1e-14;

    var maxi: usize = 0; var maxj: usize = 1; var maxv: T = 0; var phi: T = undefined; @memcpy(J.data, A.data); C.identity();

    for (0..A.rows) |i| for (i + 1..A.cols) |j| if (@abs(J.at(i, j)) > @abs(maxv)) {
        maxi = i; maxj = j; maxv = J.at(i, j);
    };

    while (@abs(maxv) > tolerance) {

        phi = 0.5 * std.math.atan(2 * maxv / (J.at(maxi, maxi) - J.at(maxj, maxj))); T1.identity();

        T1.ptr(maxi, maxi).* = std.math.cos(phi); T1.ptr(maxj, maxj).* =  T1.at(maxi, maxi);
        T1.ptr(maxj, maxi).* = std.math.sin(phi); T1.ptr(maxi, maxj).* = -T1.at(maxj, maxi);

        mat.mm(T, T2, J.*, T1.*); mat.mam(T, J, T1.*, T2.*); mat.mm(T, T2, C.*, T1.*); @memcpy(C.data, T2.data);

        maxv = 0;

        for (0..A.rows) |i| for (i + 1..A.cols) |j| if (@abs(J.at(i, j)) > @abs(maxv)) {
            maxi = i; maxj = j; maxv = J.at(i, j);
        };
    }

    for (0..A.rows) |i| for (i + 1..A.cols) |j| if (J.at(i, i) > J.at(j, j)) {

        std.mem.swap(T, J.ptr(i, i), J.ptr(j, j));

        for (0..A.rows) |k| {
            std.mem.swap(T, C.ptr(k, i), C.ptr(k, j));
        }
    };
}

/// Eigenvalue finder using QR decomposition. The eigenvalues are stored in the diagonal of the matrix J, and the eigenvectors are stored in the columns of the matrix C.
pub fn eighQr(comptime T: type, J: *Matrix(T), C: *Matrix(T), A: Matrix(T), allocator: std.mem.Allocator) !void {
    const maxiter: usize = 100000; const tolerance: T = 1e-12;

    var QT = try Matrix(T).init(A.rows, A.cols, allocator); defer QT.deinit();
    var T1 = try Matrix(T).init(A.rows, A.cols, allocator); defer T1.deinit();
    var T2 = try Matrix(T).init(A.rows, A.cols, allocator); defer T2.deinit();
    var T3 = try Vector(T).init(A.rows,         allocator); defer T3.deinit();

    try tridiagonalize(T, J, &QT, A, &T3); C.identity();

    for (0..maxiter) |i| {

        var offnorm: T = 0; for (0..A.rows) |j| for (0..A.cols) |k| {if (j != k) offnorm += J.at(j, k) * J.at(j, k);};

        if (std.math.sqrt(offnorm) < tolerance) break;

        const mu = J.at(A.rows - 1, A.cols - 1);

        for (0..J.rows) |j| J.ptr(j, j).* -= mu;

        qr(T, &T1, &T2, J.*, &T3); mat.mm(T, J, T2, T1); for (0..J.rows) |j| J.ptr(j, j).* += mu;

        mat.mm(T, &T2, C.*, T1); @memcpy(C.data, T2.data);

        if (i == maxiter - 1) return error.EighQrIterationsExceeded;
    }

    mat.mm(T, &T2, QT, C.*); @memcpy(C.data, T2.data);

    for (0..A.rows) |i| for (i + 1..A.cols) |j| if (J.at(i, i) > J.at(j, j)) {

        std.mem.swap(T, J.ptr(i, i), J.ptr(j, j));

        for (0..A.rows) |k| {
            std.mem.swap(T, C.ptr(k, i), C.ptr(k, j));
        }
    };
}

/// Solve the linear system Ax = b. The output vector is stored in the vector x. The matrix A is converted to a row-echelon form. Vector b is also modified.
pub fn linsolve(comptime T: type, x: *Vector(T), A: *Matrix(T), b: *Vector(T)) void {
    for (0..A.cols - 1) |j| for (j + 1..A.rows) |i| {

        const factor = A.at(i, j) / A.at(j, j);

        for (j..A.cols) |k| {
            A.ptr(i, k).* -= factor * A.at(j, k);
        }

        b.ptr(i).* -= factor * b.at(j);
    };

    for (0..x.rows) |i| {

        for (0..i) |j| {
            b.ptr(b.rows - i - 1).* -= A.at(A.rows - i - 1, A.cols - j - 1) * x.at(x.rows - j - 1);
        }

        x.ptr(x.rows - i - 1).* = b.at(x.rows - i - 1) / A.at(x.rows - i - 1, x.rows - i - 1);
    }
}

/// QR decomposition of a matrix A. The output matrices Q and R are stored in the matrices Q and R, respectively.
pub fn qr(comptime T: type, Q: *Matrix(T), R: *Matrix(T), A: Matrix(T), T1: *Vector(T)) void {
    const tolerance: T = 1e-12; @memcpy(R.data, A.data); Q.identity();

    for (0..A.rows) |i| {
        
        var x = T1.slice(0, A.rows - i);
        
        for (i..A.rows, 0..) |j, k| x.ptr(k).* = R.at(j, i);

        x.ptr(0).* -= (if (x.at(0) != 0) -mth.sgn(x.at(0)) else -1.0) * x.norm();

        if (x.norm() < tolerance) continue;

        vec.divs(T, &x, x, x.norm());

        for (0..A.rows - i) |j| {

            var dj: T = 0;

            for (0..A.rows - i) |k| dj += x.at(k) * R.at(i + k, i + j);

            for (0..A.rows - i) |k| R.ptr(i + k, i + j).* -= 2 * dj * x.at(k);
        }

        for (0..A.rows) |j| {

            var wj: T = 0;

            for (0..A.rows - i) |k| wj += Q.at(j, i + k) * x.at(k);

            for (0..A.rows - i) |k| Q.ptr(j, i + k).* -= 2 * wj * x.at(k);
        }
    }
}

/// Tridiagonalize the matrix A using Householder reflections into the matrix D such that A = Q * D * Q^T, where Q is an orthogonal matrix and D is a tridiagonal matrix. The output matrices Q and D are stored in the matrices Q and D, respectively.
pub fn tridiagonalize(comptime T: type, D: *Matrix(T), Q: *Matrix(T), A: Matrix(T), T1: *Vector(T)) !void {
    const tolerance: T = 1e-12; @memcpy(D.data, A.data); Q.identity(); if (A.rows < 3) return;

    for (0..A.rows - 2) |i| {

        var x = T1.slice(0, A.rows - i - 1);

        for (i + 1..A.rows, 0..) |j, k| x.ptr(k).* = D.at(j, i);

        if (x.norm() < tolerance) continue;

        x.ptr(0).* += mth.sgn(x.at(0)) * x.norm(); vec.divs(T, &x, x, x.norm());

        for (0..A.rows - i) |j| {

            var dj: T = 0;

            for (0..A.rows - i - 1) |k| dj += x.at(k) * D.at(i + k + 1, i + j);

            for (0..A.rows - i - 1) |k| D.ptr(i + k + 1, i + j).* -= 2 * dj * x.at(k);
        }

        for (0..A.rows - i) |j| {

            var dj: T = 0;

            for (0..A.rows - i - 1) |k| dj += x.at(k) * D.at(i + j, i + k + 1);

            for (0..A.rows - i - 1) |k| D.ptr(i + j, i + k + 1).* -= 2 * dj * x.at(k);
        }

        for (0..D.rows) |j| for (0..D.cols) |k| {
            if (j != k) D.ptr(j, k).* = (D.at(j, k) + D.at(k, j)) / 2;
        };

        for (0..A.rows - i) |j| {

            var wj: T = 0;

            for (0..A.rows - i - 1) |k| wj += x.at(k) * Q.at(i + j, i + k + 1);

            for (0..A.rows - i - 1) |k| Q.ptr(i + j, i + k + 1).* -= 2 * wj * x.at(k);
        }
    }
}
