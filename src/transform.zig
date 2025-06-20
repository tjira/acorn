//! Module for the integral transformation.

const std = @import("std"); const cwp = @import("cwrapper.zig");

const ten = @import("tensor.zig");

const Matrix = @import("matrix.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

/// Transforms coefficients obtained from the HF calculation from the MO basis to the MS basis.
pub fn cfsAO2AS(comptime T: type, C_MS: *Matrix(T), C_MO: Matrix(T)) void {
    C_MS.fill(0);

    for (0..C_MO.rows) |i| for (0..C_MO.cols) |j| {
        C_MS.ptr(i            , 2 * j    ).* = C_MO.at(i, j);
        C_MS.ptr(i + C_MO.rows, 2 * j + 1).* = C_MO.at(i, j);
    };
}

/// Transforms the one-electron integrals from the AO basis to the AS basis.
pub fn oneAO2AS(comptime T: type, A_AS: *Matrix(T), A_AO: Matrix(T)) void {
    A_AS.fill(0);

    for (0..A_AO.rows) |i| for (0..A_AO.cols) |j| {
        A_AS.ptr(i            , j            ).* = A_AO.at(i, j);
        A_AS.ptr(i + A_AO.rows, j + A_AO.cols).* = A_AO.at(i, j);
    };
}

/// Transforms the one-electron integrals from the AO basis to the MO basis or from the AS basis to the MS basis.
pub fn oneAO2MO(comptime T: type, A_MO: *Matrix(T), A_AO: Matrix(T), C_MO: Matrix(T)) !void {
    var A_AO_T = try Matrix(T).init(A_AO.rows, A_AO.cols, A_AO.allocator); defer A_AO_T.deinit();

    cwp.Blas(T).dgemm(&A_AO_T, A_AO, false, C_MO,   false);
    cwp.Blas(T).dgemm(A_MO,    C_MO, true,  A_AO_T, false);
}

/// Transforms the one-electron integrals from the AO basis to the MS basis.
pub fn oneAO2MS(comptime T: type, A_MS: *Matrix(T), A_AO: Matrix(T), C_MO: Matrix(T)) !void {
    var C_MS = try Matrix(T).init(2 * C_MO.rows, 2 * C_MO.cols, C_MO.allocator); defer C_MS.deinit();
    var A_AS = try Matrix(T).init(A_MS.rows, A_MS.cols,         A_MS.allocator); defer A_AS.deinit();

    cfsAO2AS(T, &C_MS, C_MO); oneAO2AS(T, A_MS, A_AO);

    cwp.Blas(T).dgemm(&A_AS, A_MS.*, false, C_MS, false);
    cwp.Blas(T).dgemm(A_MS,  C_MS,   true,  A_AS, false);
}

/// Transforms the two-electron integrals from the AO basis to the AS basis.
pub fn twoAO2AS(comptime T: type, A_AS: *Tensor(T), A_AO: Tensor(T)) void {
    A_AS.fill(0);

    for (0..A_AO.shape[0]) |i| for (0..A_AO.shape[1]) |j| for (0..A_AO.shape[2]) |k| for (0..A_AO.shape[3]) |l| {
        A_AS.ptr(&[_]usize{i,                 j,                 k                , l                }).* = A_AO.at(&[_]usize{i, j, k, l});
        A_AS.ptr(&[_]usize{i,                 j                , k + A_AO.shape[2], l + A_AO.shape[3]}).* = A_AO.at(&[_]usize{i, j, k, l});
        A_AS.ptr(&[_]usize{i + A_AO.shape[0], j + A_AO.shape[1], k,                 l                }).* = A_AO.at(&[_]usize{i, j, k, l});
        A_AS.ptr(&[_]usize{i + A_AO.shape[0], j + A_AO.shape[1], k + A_AO.shape[2], l + A_AO.shape[3]}).* = A_AO.at(&[_]usize{i, j, k, l});
    };
}

/// Transforms the two-electron integrals from the AO basis to the MO basis or from the AS basis to the MS basis.
pub fn twoAO2MO(comptime T: type, A_MO: *Tensor(T), A_AO: Tensor(T), C_MO: Matrix(T), allocator: std.mem.Allocator) !void {
    var A_AO_T = try Tensor(T).init(A_AO.shape, allocator); defer A_AO_T.deinit();

    try ten.contract(T, &A_AO_T, C_MO, &[_]usize{0}, A_AO,   &[_]usize{0}, allocator);
    try ten.contract(T, A_MO,    C_MO, &[_]usize{0}, A_AO_T, &[_]usize{1}, allocator);
    try ten.contract(T, &A_AO_T, C_MO, &[_]usize{0}, A_MO.*, &[_]usize{2}, allocator);
    try ten.contract(T, A_MO,    C_MO, &[_]usize{0}, A_AO_T, &[_]usize{3}, allocator);
}

/// Transforms the two-electron integrals from the AO basis to the MS basis.
pub fn twoAO2MS(comptime T: type, A_MS: *Tensor(T), A_AO: Tensor(T), C_MO: Matrix(T), allocator: std.mem.Allocator) !void {
    var C_MS = try Matrix(T).init(2 * C_MO.rows, 2 * C_MO.cols, allocator); defer C_MS.deinit();
    var A_AS = try Tensor(T).init(A_MS.shape,                   allocator); defer A_AS.deinit();

    cfsAO2AS(T, &C_MS, C_MO); twoAO2AS(T, A_MS, A_AO);

    try ten.contract(T, &A_AS, C_MS, &[_]usize{0}, A_MS.*, &[_]usize{0}, allocator);
    try ten.contract(T, A_MS,  C_MS, &[_]usize{0}, A_AS,   &[_]usize{1}, allocator);
    try ten.contract(T, &A_AS, C_MS, &[_]usize{0}, A_MS.*, &[_]usize{2}, allocator);
    try ten.contract(T, A_MS,  C_MS, &[_]usize{0}, A_AS,   &[_]usize{3}, allocator);
}
