const std = @import("std"); const Complex = std.math.complex.Complex;

const mat = @import("acorn").mat;

const Matrix = @import("acorn").Matrix;

test "mat_init" {
    var A: Matrix(f64) = undefined; var B: Matrix(Complex(f64)) = undefined;

    A = try @TypeOf(A).init(1, 1, std.testing.allocator); try std.testing.expect(A.rows == 1 and A.cols == 1 and A.data.len == 1); A.deinit();
    A = try @TypeOf(A).init(2, 2, std.testing.allocator); try std.testing.expect(A.rows == 2 and A.cols == 2 and A.data.len == 4); A.deinit();
    A = try @TypeOf(A).init(2, 3, std.testing.allocator); try std.testing.expect(A.rows == 2 and A.cols == 3 and A.data.len == 6); A.deinit();
    A = try @TypeOf(A).init(3, 2, std.testing.allocator); try std.testing.expect(A.rows == 3 and A.cols == 2 and A.data.len == 6); A.deinit();
    A = try @TypeOf(A).init(3, 3, std.testing.allocator); try std.testing.expect(A.rows == 3 and A.cols == 3 and A.data.len == 9); A.deinit();

    B = try @TypeOf(B).init(1, 1, std.testing.allocator); try std.testing.expect(B.rows == 1 and B.cols == 1 and B.data.len == 1); B.deinit();
    B = try @TypeOf(B).init(2, 2, std.testing.allocator); try std.testing.expect(B.rows == 2 and B.cols == 2 and B.data.len == 4); B.deinit();
    B = try @TypeOf(B).init(2, 3, std.testing.allocator); try std.testing.expect(B.rows == 2 and B.cols == 3 and B.data.len == 6); B.deinit();
    B = try @TypeOf(B).init(3, 2, std.testing.allocator); try std.testing.expect(B.rows == 3 and B.cols == 2 and B.data.len == 6); B.deinit();
    B = try @TypeOf(B).init(3, 3, std.testing.allocator); try std.testing.expect(B.rows == 3 and B.cols == 3 and B.data.len == 9); B.deinit();
}
