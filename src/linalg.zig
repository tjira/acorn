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

            if (rnorm > tolerance) {
                for (0..r.rows) |j| {r.ptr(j).* /= rnorm;} try basis.append(r);
            }
        }
    }

    for (basis.items) |v| v.deinit();
}

/// The recursive divide-and-conquer algorithm for finding the eigenvalues and eigenvectors of a real symmetric matrix A. The eigenvalues are stored in the diagonal of the matrix J, and the eigenvectors are stored in the columns of the matrix C.
pub fn eighDac(comptime T: type, J: *Matrix(T), C: *Matrix(T), A: Matrix(T), allocator: std.mem.Allocator) !void {
    const maxiter: usize = 1000000; const tol: T = 1e-14;

    var D   = try Matrix(T).init(A.rows,           A.cols,           allocator); defer   D.deinit();
    var Q   = try Matrix(T).init(A.rows,           A.cols,           allocator); defer   Q.deinit();
    var TD  = try Matrix(T).init(A.rows,           A.cols,           allocator); defer  TD.deinit();
    var TQ  = try Matrix(T).init(A.rows,           A.cols,           allocator); defer  TQ.deinit();
    var T1  = try Matrix(T).init(A.rows / 2,       A.cols / 2,       allocator); defer  T1.deinit();
    var T1J = try Matrix(T).init(A.rows / 2,       A.cols / 2,       allocator); defer T1J.deinit();
    var T1C = try Matrix(T).init(A.rows / 2,       A.cols / 2,       allocator); defer T1C.deinit();
    var T2  = try Matrix(T).init(A.rows - T1.rows, A.cols - T1.cols, allocator); defer  T2.deinit();
    var T2J = try Matrix(T).init(A.rows - T1.rows, A.cols - T1.cols, allocator); defer T2J.deinit();
    var T2C = try Matrix(T).init(A.rows - T1.rows, A.cols - T1.cols, allocator); defer T2C.deinit();
    var T3  = try Matrix(T).init(A.rows,           A.cols,           allocator); defer  T3.deinit();
    var T4  = try Matrix(T).init(A.rows,           A.cols,           allocator); defer  T4.deinit();
    var T5  = try Vector(T).init(A.rows,                             allocator); defer  T5.deinit();
    var U   = try Vector(T).init(A.rows,                             allocator); defer   U.deinit();
    var W   = try Vector(T).init(A.rows,                             allocator); defer   W.deinit();

    if (A.rows == 1) {J.ptr(0, 0).* = A.at(0, 0); C.ptr(0, 0).* = 1; return;}

    try tridiagonalize(T, &TD, &TQ, A, &T5); C.fill(0);

    for (0..T1.rows) |i| for (0..T1.cols) |j| {T1.ptr(i, j).* = TD.at(i,           j          );};
    for (0..T2.rows) |i| for (0..T2.cols) |j| {T2.ptr(i, j).* = TD.at(i + T1.rows, j + T1.cols);};

    const beta = TD.at(A.rows / 2 - 1, A.cols / 2); T1.ptr(T1.rows - 1, T1.cols - 1).* -= beta; T2.ptr(0, 0).* -= beta;

    try eighDac(T, &T1J, &T1C, T1, allocator);
    try eighDac(T, &T2J, &T2C, T2, allocator);

    for (0..T1.rows) |i| for (0..T1.cols) |j| {D.ptr(i,           j          ).* = T1J.at(i, j);};
    for (0..T2.rows) |i| for (0..T2.cols) |j| {D.ptr(i + T1.rows, j + T1.cols).* = T2J.at(i, j);};
    for (0..T1.rows) |i| for (0..T1.cols) |j| {Q.ptr(i,           j          ).* = T1C.at(i, j);};
    for (0..T2.rows) |i| for (0..T2.cols) |j| {Q.ptr(i + T1.rows, j + T1.cols).* = T2C.at(i, j);};

    for (0..T1.cols) |i| W.ptr(i          ).* = T1C.at(T1.rows - 1, i);
    for (0..T2.cols) |i| W.ptr(i + T1.rows).* = T2C.at(0,           i);

    for (0..A.rows) |i| for (i + 1..A.cols) |j| if (D.at(i, i) > D.at(j, j)) {

        std.mem.swap(T, D.ptr(i, i), D.ptr(j, j));
        std.mem.swap(T, W.ptr(i   ), W.ptr(j   ));

        for (0..A.rows) |k| {
            std.mem.swap(T, Q.ptr(k, i), Q.ptr(k, j));
        }
    };

    const SecularD0 = struct { fn get (lambda: T, rho: T, WW: Vector(T), DD: Matrix(T)) T {
        var sum: T = 0; for (0..WW.rows) |i| sum += WW.at(i) * WW.at(i) / (DD.at(i, i) - lambda); return 1 + rho * sum;
    }};
    const SecularD1 = struct { fn get (lambda: T, rho: T, WW: Vector(T), DD: Matrix(T)) T {
        var sum: T = 0; for (0..WW.rows) |i| sum += WW.at(i) * WW.at(i) / ((DD.at(i, i) - lambda) * (DD.at(i, i) - lambda)); return rho * sum;
    }};

    for (0..A.rows) |i| {

        var a: T = 0; var b: T = 0;

        if (beta < 0 and i == 0         ) {a = D.at(0,          0         ) - 1e6; b = D.at(0,          0         )      ;}
        if (beta < 0 and i != 0         ) {a = D.at(i - 1,      i - 1     )      ; b = D.at(i,          i         )      ;}
        if (beta > 0 and i == A.rows - 1) {a = D.at(A.rows - 1, A.rows - 1)      ; b = D.at(A.rows - 1, A.rows - 1) + 1e6;}
        if (beta > 0 and i != A.rows - 1) {a = D.at(i,          i         )      ; b = D.at(i + 1,      i + 1     )      ;}

        var x = (a + b) / 2;

        if (@abs(b - a) > tol) for (0..maxiter) |j| {

            const fx  = SecularD0.get(x, beta, W, D);
            const dfx = SecularD1.get(x, beta, W, D);

            if (fx * beta > 0) {b = x;} else a = x; const xp = x;

            if (@abs(dfx) > tol) {x = x - fx / dfx;} else x = (a + b) / 2;

            if (x < a or x > b) {x = (a + b) / 2;}

            if (@abs(x - xp) < tol) break;

            if (j == maxiter - 1) return error.EighDacSecularIterationsExceeded;
        };

        J.ptr(i, i).* = x;
    }

    for (0..A.rows) |i| {

        for (0..U.rows) |j| U.ptr(j).* = W.at(j) / (D.at(j, j) - J.at(i, i) + 1e-14);

        vec.divs(T, &U, U, U.norm());

        for (0..A.rows) |j| {
            for (0..A.rows) |k| {
                C.ptr(j, i).* += Q.at(j, k) * U.at(k);
            }
        }
    }

    mat.mm(T, &T3, TQ, C.*); @memcpy(C.data, T3.data);
}

/// Find the eigenvalues and eigenvectors of a real symmetric matrix A. The eigenvalues are stored in the diagonal of the matrix J, and the eigenvectors are stored in the columns of the matrix C. The matrices T1 and T2 are temporary matrices used in the computation.
pub fn eighJacobi(comptime T: type, J: *Matrix(T), C: *Matrix(T), A: Matrix(T), T1: *Matrix(T), T2: *Matrix(T)) void {
    const tol: T = 1e-14;

    var maxi: usize = 0; var maxj: usize = 1; var maxv: T = 0; var phi: T = undefined; @memcpy(J.data, A.data); C.identity();

    for (0..A.rows) |i| for (i + 1..A.cols) |j| if (@abs(J.at(i, j)) > @abs(maxv)) {
        maxi = i; maxj = j; maxv = J.at(i, j);
    };

    while (@abs(maxv) > tol) {

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

        for (0..A.rows) |j| {

            var wj: T = 0;

            for (0..A.rows - i - 1) |k| wj += x.at(k) * Q.at(j, i + k + 1);

            for (0..A.rows - i - 1) |k| Q.ptr(j, i + k + 1).* -= 2 * wj * x.at(k);
        }
    }
}
