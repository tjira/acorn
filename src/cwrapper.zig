//! Wrappers for the C/C++ libraries.

const cblas   = @cImport(@cInclude("cblas.h"  ));
const eigen   = @cImport(@cInclude("eigen.h"  ));
const exprtk  = @cImport(@cInclude("exprtk.h" ));
const lapacke = @cImport(@cInclude("lapacke.h"));
const libint  = @cImport(@cInclude("libint.h" ));

const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const System = @import("system.zig").System;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("vector.zig").Vector;

/// Blas struct.
pub fn Blas(comptime T: type) type {
    return struct {

        /// Calculate the matrix-matrix product C = A * B and store the result in C. The AT and BT flags control whether the matrices A and B are transposed or not.
        pub fn dgemm(C: *Matrix(T), A: Matrix(T), AT: bool, B: Matrix(T), BT: bool) !void {
            if (!AT and !BT and (C.rows != A.rows or C.cols != B.cols or A.cols != B.rows)) return error.ErrorInMatrixDimensions;
            if ( AT and !BT and (C.rows != A.cols or C.cols != B.cols or A.rows != B.rows)) return error.ErrorInMatrixDimensions;
            if (!AT and  BT and (C.rows != A.rows or C.cols != B.rows or A.cols != B.cols)) return error.ErrorInMatrixDimensions;
            if ( AT and  BT and (C.rows != A.cols or C.cols != B.rows or A.rows != B.cols)) return error.ErrorInMatrixDimensions;

            const m: i32 = @intCast(C.rows); const n: i32 = @intCast(C.cols); const k: i32 = @intCast(if (!AT) A.cols else A.rows);

            const at: c_uint = if (AT) cblas.CblasTrans else cblas.CblasNoTrans;
            const bt: c_uint = if (BT) cblas.CblasTrans else cblas.CblasNoTrans;

            const LDA = if (!AT) k else m;
            const LDB = if (!BT) n else k;

            cblas.cblas_dgemm(cblas.CblasRowMajor, at, bt, m, n, k, 1.0, &A.data[0], LDA, &B.data[0], LDB, 0.0, &C.data[0], n);
        }

        /// Calculate the matrix-matrix product C = A * B and store the result in C. The AH and BH flags control whether the matrices A and B are Hermitian transposed or not.
        pub fn zgemm(C: *Matrix(std.math.Complex(T)), A: Matrix(std.math.Complex(T)), AH: bool, B: Matrix(std.math.Complex(T)), BH: bool) !void {
            if (C.rows != A.rows or C.cols != B.cols or A.cols != B.rows) return error.ErrorInMatrixDimensions;

            const m: i32 = @intCast(C.rows); const n: i32 = @intCast(C.cols); const k: i32 = @intCast(A.cols);

            const ah: c_uint = if (AH) cblas.CblasConjTrans else cblas.CblasNoTrans;
            const bh: c_uint = if (BH) cblas.CblasConjTrans else cblas.CblasNoTrans;

            const LDA = if (!AH) k else m;
            const LDB = if (!BH) n else k;

            cblas.cblas_zgemm(cblas.CblasRowMajor, ah, bh, m, n, k, &std.math.Complex(T).init(1.0, 0.0), &A.data[0], LDA, &B.data[0], LDB, &std.math.Complex(T).init(0.0, 0.0), &C.data[0], n);
        }
    };
}

/// Eigen struct.
pub fn Eigen(comptime T: type) type {
    return struct {

        /// Determinant of a matrix.
        pub fn det(A: Matrix(T)) T {
            return eigen.det(&A.data[0], A.rows);
        }

        /// Logarithm of a matrix.
        pub fn logm(B: *Matrix(T), A: Matrix(T)) void {
            eigen.logm(&B.data[0], &A.data[0], A.rows);
        }
    };
}

/// Expression struct.
pub fn Expression(comptime T: type) type {
    return struct {
        expression: *anyopaque,

        /// Initializer function.
        pub fn init(string: []const u8, nvars: usize, allocator: std.mem.Allocator) !Expression(T) {
            const buffer = try allocator.dupeZ(u8, string); defer allocator.free(buffer); return .{.expression = exprtk.compile(buffer, nvars).?};
        }

        /// Deinitialize the expression.
        pub fn deinit(self: Expression(T)) void {
            exprtk.deinit(self.expression);
        }

        /// Evaluate the compiled expression at the specified point.
        pub fn evaluate(self: Expression(T), vars: Vector(T), time: T) T {
            return exprtk.evaluate(self.expression, vars.data.ptr, time);
        }
    };
}

/// Lapack struct.
pub fn Lapack(comptime T: type) type {
    return struct {

        /// Calculate the condition number of a matrix A. The result returned.
        pub fn dgecon(A: Matrix(T)) !T {
            const n: i32 = @intCast(A.rows); var rcond: T = 0;

            var AP = try Matrix(T).init(@as(usize, @intCast(n)), @as(usize, @intCast(n)), A.allocator); var ALU = try A.clone();

            try dgetrf(&AP, &ALU, A);

            const info = lapacke.LAPACKE_dgecon(lapacke.LAPACK_ROW_MAJOR, '1', n, &ALU.data[0], n, A.onorm(), &rcond);

            if (info != 0) return error.ErrorInConditionNumberCalculation;

            AP.deinit(); ALU.deinit();

            return rcond;
        }

        /// Calculate the Schur decomposition of matrix A.
        pub fn dgees(D: *Matrix(T), Q: *Matrix(T), A: Matrix(T)) !void {
            const n: i32 = @intCast(A.rows); A.memcpy(D.*);

            var t1: i32 = undefined; var T1 = try Vector(T).init(@as(usize, @intCast(n)), A.allocator); var T2 = try Vector(T).init(@as(usize, @intCast(n)), A.allocator);

            const info = lapacke.LAPACKE_dgees(lapacke.LAPACK_ROW_MAJOR, 'V', 'N', null, n, &D.data[0], n, &t1, &T1.data[0], &T2.data[0], &Q.data[0], n);

            if (info != 0) return error.ErrorInRealSchurDecomposition;

            T1.deinit(); T2.deinit();
        }

        /// Perform the SVD decomposition of a general matrix.
        pub fn dgesdd(U: *Matrix(T), S: *Matrix(T), VT: *Matrix(T), A: Matrix(T)) !void {
            const m: i32 = @intCast(A.rows); const n: i32 = @intCast(A.cols);

            var T1 = try A.clone();

            const info = lapacke.LAPACKE_dgesdd(lapacke.LAPACK_ROW_MAJOR, 'A', m, n, &T1.data[0], n, &S.data[0], &U.data[0], m, &VT.data[0], n);

            for (0..S.cols) |i| std.mem.swap(T, S.ptr(i, i), S.ptr(0, i));

            if (info != 0) return error.ErrorInSingularValueDecomposition;

            T1.deinit();
        }

        /// Solve the linear system Ax = b using the LU decomposition of A. The result is stored in x.
        pub fn dgesv(x: *Vector(T), A: Matrix(T), b: Vector(T)) !void {
            const n: i32 = @intCast(A.rows); b.memcpy(x.*);

            var T1 = try A.clone(); var T2 = try Vector(i32).init(@as(usize, @intCast(n)), A.allocator);

            const info = lapacke.LAPACKE_dgesv(lapacke.LAPACK_ROW_MAJOR, n, 1, &T1.data[0], n, &T2.data[0], &x.data[0], 1);

            if (info != 0) return error.ErrorInLinearSystemSolution;

            T1.deinit(); T2.deinit();
        }

        /// Perform the LU factorization of a general matrix A. The permutation matrix P is stored in the diagonal of P and the LU decomposition is stored in LU.
        pub fn dgetrf(P: *Matrix(T), LU: *Matrix(T), A: Matrix(T)) !void {
            const m: i32 = @intCast(A.rows); const n: i32 = @intCast(A.cols); A.memcpy(LU.*); P.fill(0);

            var T1 = try Vector(i32).init(@as(usize, @intCast(@min(m, n))), A.allocator);

            const info = lapacke.LAPACKE_dgetrf(lapacke.LAPACK_ROW_MAJOR, m, n, &LU.data[0], n, &T1.data[0]);

            for (0..T1.rows) |i| P.ptr(i, @as(usize, @intCast(T1.at(i) - 1))).* = 1;

            if (info != 0) return error.ErrorInLUFactorization;

            T1.deinit();
        }

        /// Calculate the eigenvalues and eigenvectors of a symmetric matrix A. The eigenvalues are stored in J and the eigenvectors are stored in C.
        pub fn dsyevd(J: *Matrix(T), C: *Matrix(T), A: Matrix(T)) !void {
            const n: i32 = @intCast(A.rows); A.memcpy(C.*); J.fill(0);

            const info = lapacke.LAPACKE_dsyevd(lapacke.LAPACK_ROW_MAJOR, 'V', 'U', n, &C.data[0], n, &J.data[0]);

            for (0..J.rows) |i| std.mem.swap(T, J.ptr(i, i), J.ptr(0, i));

            if (info != 0) return error.ErrorInEigenvalueCalculation;
        }

        /// Finds the eigenvalues and eigenvectors of a symmetric-definite generalized eigenproblem A*x = λ*B*x. The eigenvalues are stored in J and the eigenvectors are stored in C.
        pub fn dsygvd(J: *Matrix(T), C: *Matrix(T), A: Matrix(T), B: Matrix(T)) !void {
            const n: i32 = @intCast(A.rows); A.memcpy(C.*); J.fill(0);

            var T1 = try B.clone();

            const info = lapacke.LAPACKE_dsygvd(lapacke.LAPACK_ROW_MAJOR, 1, 'V', 'U', n, &C.data[0], n, &T1.data[0], n, &J.data[0]);

            for (0..J.rows) |i| std.mem.swap(T, J.ptr(i, i), J.ptr(0, i));

            if (info != 0) return error.ErrorInEigenvalueCalculation;

            T1.deinit();
        }
    };
}

/// Libint struct.
pub fn Libint(comptime T: type) type {
    return struct {

        /// Calculate the Coulomb integrals.
        pub fn coulomb(I: *Tensor(T), system: System(T), basis: std.ArrayList(T)) void {
            libint.coulomb(&I.data[0], system.atoms.rows, &system.atoms.data[0], &system.coords.data[0], basis.items.len, basis.items.ptr);
        }

        /// Calculate the kinetic integrals.
        pub fn kinetic(I: *Matrix(T), system: System(T), basis: std.ArrayList(T)) void {
            libint.kinetic(&I.data[0], system.atoms.rows, &system.atoms.data[0], &system.coords.data[0], basis.items.len, basis.items.ptr);
        }

        /// Calculate the nuclear integrals.
        pub fn nuclear(I: *Matrix(T), system: System(T), basis: std.ArrayList(T)) void {
            libint.nuclear(&I.data[0], system.atoms.rows, &system.atoms.data[0], &system.coords.data[0], basis.items.len, basis.items.ptr);
        }

        /// Calculate the overlap integrals.
        pub fn overlap(I: *Matrix(T), system: System(T), basis: std.ArrayList(T)) void {
            libint.overlap(&I.data[0], system.atoms.rows, &system.atoms.data[0], &system.coords.data[0], basis.items.len, basis.items.ptr);
        }

        /// Calculate the Coulomb integrals and contract them with density.
        pub fn fock(F: *Matrix(T), system: System(T), basis: std.ArrayList(T), D: Matrix(T), generalized: bool) void {
            libint.fock(&F.data[0], system.atoms.rows, &system.atoms.data[0], &system.coords.data[0], basis.items.len, basis.items.ptr, &D.data[0], generalized);
        }
    };
}
