const std = @import("std"); const Complex = std.math.complex.Complex;

const mat = @import("acorn").mat;
const vec = @import("acorn").vec;

const Matrix = @import("acorn").Matrix;
const Vector = @import("acorn").Vector;

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

test "mat_linsolve" {
    var A: Matrix(f64) = undefined;
    var b: Vector(f64) = undefined;
    var x: Vector(f64) = undefined;
    var y: Vector(f64) = undefined;

    A = try @TypeOf(A).init(2, 2, std.testing.allocator);
    b = try @TypeOf(b).init(2,    std.testing.allocator);
    x = try @TypeOf(x).init(2,    std.testing.allocator);
    y = try @TypeOf(x).init(2,    std.testing.allocator);

    A.ptr(0, 0).* = 5; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 1; A.ptr(1, 1).* = 4;

    b.ptr(0).* = 3; b.ptr(1).* = 2;

    y.ptr(0).* = 4.0 / 9.0; y.ptr(1).* = 7.0 / 18.0;

    mat.linsolve(f64, &x, A, b);

    try std.testing.expect(vec.eq(@TypeOf(x.at(0)), x, y, 1e-12));

    A.deinit(); b.deinit(); x.deinit(); y.deinit();

    A = try @TypeOf(A).init(3, 3, std.testing.allocator);
    b = try @TypeOf(b).init(3,    std.testing.allocator);
    x = try @TypeOf(x).init(3,    std.testing.allocator);
    y = try @TypeOf(x).init(3,    std.testing.allocator);

    A.ptr(0, 0).* = 25; A.ptr(0, 1).* = 2; A.ptr(0, 2).* = 1;
    A.ptr(1, 0).* = 1; A.ptr(1, 1).* = 31; A.ptr(1, 2).* = 2;
    A.ptr(2, 0).* = 2; A.ptr(2, 1).* = 3; A.ptr(2, 2).* = 18;

    b.ptr(0).* = 3; b.ptr(1).* = 2; b.ptr(2).* = 1;

    y.ptr(0).* = 521.0 / 4571.0; y.ptr(1).* = 115.0 / 1959.0; y.ptr(2).* = 454.0 / 13713.0;

    mat.linsolve(f64, &x, A, b);

    try std.testing.expect(vec.eq(@TypeOf(x.at(0)), x, y, 1e-12));

    A.deinit(); b.deinit(); x.deinit(); y.deinit();
}
