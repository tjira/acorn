const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;

const allocator = std.heap.page_allocator; const tolerance: f64 = 1e-12;

test "matrixAdd" {
    const A = try Matrix(f64).init(2, 2, &[_]f64{1, 2, 3, 4}, allocator); // define the first matrix
    const B = try Matrix(f64).init(2, 2, &[_]f64{2, 3, 4, 5}, allocator); // define the second matrix
    const result = try A.add(B); // calculate their sum
    const expected = try Matrix(f64).init(2, 2, &[_]f64{3, 5, 7, 9}, allocator); // define the expected result
    try std.testing.expect(try result.eqApprox(expected, tolerance)); // return the result of the test
}
