const std = @import("std"); const Complex = std.math.complex.Complex;

const mat = @import("acorn").mat;

const Matrix = @import("acorn").Matrix;

const sqrt1_2 = std.math.sqrt1_2;

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

test "mat_eigh" {
    var A:  Matrix(f64) = undefined;
    var B:  Matrix(f64) = undefined;
    var C:  Matrix(f64) = undefined;
    var AJ: Matrix(f64) = undefined;
    var AC: Matrix(f64) = undefined;
    var T1: Matrix(f64) = undefined;
    var T2: Matrix(f64) = undefined;

    A  = try @TypeOf(A ).init(2, 2, std.testing.allocator);
    B  = try @TypeOf(B ).init(2, 2, std.testing.allocator);
    C  = try @TypeOf(C ).init(2, 2, std.testing.allocator);
    AJ = try @TypeOf(AJ).init(2, 2, std.testing.allocator);
    AC = try @TypeOf(AC).init(2, 2, std.testing.allocator);
    T1 = try @TypeOf(T1).init(2, 2, std.testing.allocator);
    T2 = try @TypeOf(T2).init(2, 2, std.testing.allocator);

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 2; A.ptr(1, 1).* = 1;

    B.ptr(0, 0).* = -1; B.ptr(0, 1).* = 0;
    B.ptr(1, 0).* =  0; B.ptr(1, 1).* = 3;

    C.ptr(0, 0).* = -sqrt1_2; C.ptr(0, 1).* = sqrt1_2;
    C.ptr(1, 0).* =  sqrt1_2; C.ptr(1, 1).* = sqrt1_2;

    mat.eigh(@TypeOf(A.at(0, 0)), &AJ, &AC, A, &T1, &T2);

    try std.testing.expect(mat.eq(@TypeOf(A.at(0, 0)), AJ, B, 1e-12) and mat.eq(@TypeOf(A.at(0, 0)), AC, C, 1e-12));

    A.deinit(); B.deinit(); C.deinit(); AJ.deinit(); AC.deinit(); T1.deinit(); T2.deinit();
}
