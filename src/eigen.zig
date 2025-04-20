//! Wrapper for the Eigen C++ library.

const eigen = @cImport(@cInclude("eigen.h"));

const std = @import("std");

const Tensor = @import("tensor.zig").Tensor;

pub fn contract(C: anytype, A: anytype, B: anytype, pairs: []const i32) !void {
    C.fill(0); var dC: []const usize = undefined; var dA: []const usize = undefined; var dB: []const usize = undefined;

    if (comptime std.meta.fieldIndex(@TypeOf(C.*), "shape") != null) {dC = C.shape;} else {dC = &[_]usize{C.rows, C.cols};}
    if (comptime std.meta.fieldIndex(@TypeOf(A  ), "shape") != null) {dA = A.shape;} else {dA = &[_]usize{A.rows, A.cols};}
    if (comptime std.meta.fieldIndex(@TypeOf(B  ), "shape") != null) {dB = B.shape;} else {dB = &[_]usize{B.rows, B.cols};}

    if (dA.len == 2 and dB.len == 2 and pairs.len == 2) {return eigen.contract_221(&C.data[0], &dC[0], &A.data[0], &dA[0], &B.data[0], &dB[0], &pairs[0]);}
    if (dA.len == 2 and dB.len == 4 and pairs.len == 2) {return eigen.contract_241(&C.data[0], &dC[0], &A.data[0], &dA[0], &B.data[0], &dB[0], &pairs[0]);}
    if (dA.len == 2 and dB.len == 4 and pairs.len == 4) {return eigen.contract_242(&C.data[0], &dC[0], &A.data[0], &dA[0], &B.data[0], &dB[0], &pairs[0]);}

    else return error.ContractionNotSupported;
}
