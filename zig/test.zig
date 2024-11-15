const std = @import("std");

const cdn = @import("classicaldynamics.zig");
const mpt = @import("modelpotential.zig"   );
const mat = @import("matrix.zig"           );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const allocator = std.testing.allocator; const tol = 1e-14;

test "matrixFill2x2" {
    var A = try Matrix(f64).init(2, 2, allocator); defer A.deinit(); 
    var O = try Matrix(f64).init(2, 2, allocator); defer O.deinit();
    A.fill(1); O.set(&[_]f64{1, 1, 1, 1});
    try std.testing.expect(A.eq(O, tol));
}

test "matrixHjoin2x2-2x2" {
    var A = try Matrix(f64).init(2, 2, allocator); defer A.deinit(); 
    var B = try Matrix(f64).init(2, 2, allocator); defer B.deinit(); 
    var C = try Matrix(f64).init(2, 4, allocator); defer C.deinit(); 
    var O = try Matrix(f64).init(2, 4, allocator); defer O.deinit();
    A.set(&[_]f64{1, 2, 3, 4}); B.set(&[_]f64{5, 6, 7, 8}); A.hjoin(&C, B); O.set(&[_]f64{1, 2, 5, 6, 3, 4, 7, 8});
    try std.testing.expect(C.eq(O, tol));
}

test "matrixLinspace2x2" {
    var A = try Matrix(f64).init(2, 2, allocator); defer A.deinit(); 
    var O = try Matrix(f64).init(2, 2, allocator); defer O.deinit();
    A.linspace(-2, 1); O.set(&[_]f64{-2, -1, 0, 1});
    try std.testing.expect(A.eq(O, tol));
}

test "matrixEigh2x2" {
    var A  = try Matrix(f64).init(2, 2, allocator); defer  A.deinit();
    var B  = try Matrix(f64).init(2, 2, allocator); defer  B.deinit();
    var C  = try Matrix(f64).init(2, 2, allocator); defer  C.deinit();
    var O1 = try Matrix(f64).init(2, 2, allocator); defer O1.deinit();
    var O2 = try Matrix(f64).init(2, 2, allocator); defer O2.deinit();
    A.set(&[_]f64{1, 2, 2, 1}); try mat.eigh(f64, &B, &C, A, tol); O1.set(&[_]f64{-1, 0, 0, 3}); O2.set(&[_]f64{-std.math.sqrt1_2, std.math.sqrt1_2, std.math.sqrt1_2, std.math.sqrt1_2});
    try std.testing.expect(B.eq(O1, tol) and C.eq(O2, tol));
}

test "matrixMM2x2-2x2" {
    var A = try Matrix(f64).init(2, 2, allocator); defer A.deinit();
    var B = try Matrix(f64).init(2, 2, allocator); defer B.deinit();
    var C = try Matrix(f64).init(2, 2, allocator); defer C.deinit();
    var O = try Matrix(f64).init(2, 2, allocator); defer O.deinit();
    A.set(&[_]f64{1, 2, 3, 4}); B.set(&[_]f64{5, 6, 7, 8}); mat.mm(f64, &C, A, B); O.set(&[_]f64{19, 22, 43, 50});
    try std.testing.expect(C.eq(O, tol));
}

test "matrixTranspose2x2" {
    var A = try Matrix(f64).init(2, 2, allocator); defer A.deinit();
    var B = try Matrix(f64).init(2, 2, allocator); defer B.deinit();
    var O = try Matrix(f64).init(2, 2, allocator); defer O.deinit();
    A.set(&[_]f64{1, 2, 3, 4}); mat.transpose(f64, &B, A); O.set(&[_]f64{1, 3, 2, 4});
    try std.testing.expect(B.eq(O, tol));
}
