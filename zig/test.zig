const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;

const allocator = std.testing.allocator; const tolerance: f64 = 1e-12;

const matrix_2x2_1 = [_]f64{1, 2, 5, 4};
const matrix_2x2_2 = [_]f64{1, 2, 1, 2};

test "matrix_add" {
    const A = try Matrix(f64).init(2, 2, &matrix_2x2_1, allocator); defer A.deinit();
    const B = try Matrix(f64).init(2, 2, &matrix_2x2_2, allocator); defer B.deinit();
    const result = try A.add(B); defer result.deinit();
    const expected = try Matrix(f64).init(2, 2, &[_]f64{2, 4, 6, 6}, allocator); defer expected.deinit();
    try std.testing.expect(result.eqApprox(expected, tolerance));
}

test "matrix_div" {
    const A = try Matrix(f64).init(2, 2, &matrix_2x2_1, allocator); defer A.deinit();
    const B = try Matrix(f64).init(2, 2, &matrix_2x2_2, allocator); defer B.deinit();
    const result = try A.div(B); defer result.deinit();
    const expected = try Matrix(f64).init(2, 2, &[_]f64{1, 1, 5, 2}, allocator); defer expected.deinit();
    try std.testing.expect(result.eqApprox(expected, tolerance));
}

test "matrix_mul" {
    const A = try Matrix(f64).init(2, 2, &matrix_2x2_1, allocator); defer A.deinit();
    const B = try Matrix(f64).init(2, 2, &matrix_2x2_2, allocator); defer B.deinit();
    const result = try A.mul(B); defer result.deinit();
    const expected = try Matrix(f64).init(2, 2, &[_]f64{1, 4, 5, 8}, allocator); defer expected.deinit();
    try std.testing.expect(result.eqApprox(expected, tolerance));
}

test "matrix_sub" {
    const A = try Matrix(f64).init(2, 2, &matrix_2x2_1, allocator); defer A.deinit();
    const B = try Matrix(f64).init(2, 2, &matrix_2x2_2, allocator); defer B.deinit();
    const result = try A.sub(B); defer result.deinit();
    const expected = try Matrix(f64).init(2, 2, &[_]f64{0, 0, 4, 2}, allocator); defer expected.deinit();
    try std.testing.expect(result.eqApprox(expected, tolerance));
}
