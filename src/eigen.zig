//! Wrapper for the Eigen C++ library.

const eigen = @cImport(@cInclude("eigen.h"));

const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

pub fn contract(C: anytype, A: anytype, B: anytype, pairs: []const i32) !void {
    C.fill(0);

    var dC: []const usize = undefined; var rC: usize = undefined;
    var dA: []const usize = undefined; var rA: usize = undefined;
    var dB: []const usize = undefined; var rB: usize = undefined;

    if (comptime std.meta.fieldIndex(@TypeOf(C.*), "shape") != null) {dC = C.shape; rC = C.shape.len;} else {dC = &[_]usize{C.rows, C.cols}; rC = 2;}
    if (comptime std.meta.fieldIndex(@TypeOf(A  ), "shape") != null) {dA = A.shape; rA = A.shape.len;} else {dA = &[_]usize{A.rows, A.cols}; rA = 2;}
    if (comptime std.meta.fieldIndex(@TypeOf(B  ), "shape") != null) {dB = B.shape; rB = B.shape.len;} else {dB = &[_]usize{B.rows, B.cols}; rB = 2;}

    return eigen.contract(&C.data[0], &dC[0], rC, &A.data[0], &dA[0], rA, &B.data[0], &dB[0], rB, &pairs[0], pairs.len / 2);
}

pub fn logm(B: *Matrix(f64), A: Matrix(f64)) void {
    eigen.logm(&B.data[0], &A.data[0], A.rows);
}
