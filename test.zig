const std = @import("std"); const Complex = std.math.Complex;

const cdn = @import("src/classicaldynamics.zig");
const ftr = @import("src/fouriertransform.zig" );
const mat = @import("src/matrix.zig"           );
const mpt = @import("src/modelpotential.zig"   );
const qdn = @import("src/quantumdynamics.zig"  );
const vec = @import("src/vector.zig"           );

const Vector = @import("src/vector.zig").Vector;
const Matrix = @import("src/matrix.zig").Matrix;

const CDO = @import("src/classicaldynamics.zig").ClassicalDynamicsOptions;
const QDO = @import("src/quantumdynamics.zig"  ).QuantumDynamicsOptions  ;
const MPO = @import("src/modelpotential.zig"   ).ModelPotentialOptions   ;

// VECTOR ==============================================================================================================================================================================================

test "vector_init_1" {
    var u = try Vector(f64).init(1, std.testing.allocator); defer u.deinit();

    try std.testing.expect(u.data.len == 1 and u.rows == 1);
}

test "vector_init_2" {
    var u = try Vector(f64).init(2, std.testing.allocator); defer u.deinit();

    try std.testing.expect(u.data.len == 2 and u.rows == 2);
}

test "vector_clone_1" {
    var u = try Vector(f64).init(1, std.testing.allocator); defer u.deinit();

    u.ptr(0).* = 3;

    var v = try u.clone(); defer v.deinit();

    try std.testing.expect(vec.eq(f64, u, v, 1e-12));
}

test "vector_clone_2" {
    var u = try Vector(f64).init(2, std.testing.allocator); defer u.deinit();

    u.ptr(0).* = 3;
    u.ptr(1).* = 4;

    var v = try u.clone(); defer v.deinit();

    try std.testing.expect(vec.eq(f64, u, v, 1e-12));
}

test "vector_complex_1" {
    var u = try Vector(f64         ).init(1, std.testing.allocator); defer u.deinit();
    var r = try Vector(Complex(f64)).init(1, std.testing.allocator); defer r.deinit();

    u.ptr(0).* = 3;

    r.ptr(0).* = Complex(f64).init(3, 0);

    var v = try u.complex(); defer v.deinit();

    try std.testing.expect(vec.ceq(f64, v, r, 1e-12));
}

test "vector_complex_2" {
    var u = try Vector(f64         ).init(2, std.testing.allocator); defer u.deinit();
    var r = try Vector(Complex(f64)).init(2, std.testing.allocator); defer r.deinit();

    u.ptr(0).* = 3;
    u.ptr(1).* = 4;

    r.ptr(0).* = Complex(f64).init(3, 0);
    r.ptr(1).* = Complex(f64).init(4, 0);

    var v = try u.complex(); defer v.deinit();

    try std.testing.expect(vec.ceq(f64, v, r, 1e-12));
}

test "vector_at_ptr_1" {
    var u = try Vector(f64).init(1, std.testing.allocator); defer u.deinit();

    u.ptr(0).* = 3;

    try std.testing.expect(u.at(0) == 3);
}

test "vector_at_ptr_2" {
    var u = try Vector(f64).init(2, std.testing.allocator); defer u.deinit();

    u.ptr(0).* = 3;
    u.ptr(1).* = 4;

    try std.testing.expect(u.at(0) == 3 and u.at(1) == 4);
}

test "vector_matrixptr_1" {
    var u = try Vector(f64).init(1, std.testing.allocator); defer u.deinit();

    u.ptr(0).* = 3;

    const v = u.matrixptr().vectorptr();

    try std.testing.expect(vec.eq(f64, u, v, 1e-12));
}

test "vector_matrixptr_2" {
    var u = try Vector(f64).init(2, std.testing.allocator); defer u.deinit();

    u.ptr(0).* = 3;
    u.ptr(1).* = 4;

    const v = u.matrixptr().vectorptr();

    try std.testing.expect(vec.eq(f64, u, v, 1e-12));
}

test "vector_fill_1" {
    var u = try Vector(f64).init(1, std.testing.allocator); defer u.deinit();
    var r = try Vector(f64).init(1, std.testing.allocator); defer r.deinit();

    r.ptr(0).* = 3;

    u.fill(3);

    try std.testing.expect(vec.eq(f64, u, r, 1e-12));
}

test "vector_fill_2" {
    var u = try Vector(f64).init(2, std.testing.allocator); defer u.deinit();
    var r = try Vector(f64).init(2, std.testing.allocator); defer r.deinit();

    r.ptr(0).* = 3;
    r.ptr(1).* = 3;

    u.fill(3);

    try std.testing.expect(vec.eq(f64, u, r, 1e-12));
}

test "eq_1" {
    var u = try Vector(f64).init(1, std.testing.allocator); defer u.deinit();
    var v = try Vector(f64).init(1, std.testing.allocator); defer v.deinit();

    u.ptr(0).* = 3;

    v.ptr(0).* = 3;

    try std.testing.expect(vec.eq(f64, u, v, 1e-12));
}

test "eq_2" {
    var u = try Vector(f64).init(2, std.testing.allocator); defer u.deinit();
    var v = try Vector(f64).init(2, std.testing.allocator); defer v.deinit();

    u.ptr(0).* = 3;
    u.ptr(1).* = 4;

    v.ptr(0).* = 3;
    v.ptr(1).* = 4;

    try std.testing.expect(vec.eq(f64, u, v, 1e-12));
}

test "ceq_1" {
    var u = try Vector(Complex(f64)).init(1, std.testing.allocator); defer u.deinit();
    var v = try Vector(Complex(f64)).init(1, std.testing.allocator); defer v.deinit();

    u.ptr(0).* = Complex(f64).init(3, 1);

    v.ptr(0).* = Complex(f64).init(3, 1);

    try std.testing.expect(vec.ceq(f64, u, v, 1e-12));
}

test "ceq_2" {
    var u = try Vector(Complex(f64)).init(2, std.testing.allocator); defer u.deinit();
    var v = try Vector(Complex(f64)).init(2, std.testing.allocator); defer v.deinit();

    u.ptr(0).* = Complex(f64).init(3, 1);
    u.ptr(1).* = Complex(f64).init(4, 2);

    v.ptr(0).* = Complex(f64).init(3, 1);
    v.ptr(1).* = Complex(f64).init(4, 2);

    try std.testing.expect(vec.ceq(f64, u, v, 1e-12));
}

// MATRIX ==============================================================================================================================================================================================

test "matrix_init_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();

    try std.testing.expect(A.rows == 1 and A.cols == 1 and A.data.len == 1);
}

test "matrix_init_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();

    try std.testing.expect(A.rows == 2 and A.cols == 2 and A.data.len == 4);
}

test "matrix_clone_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();

    var B = try A.clone(); defer B.deinit();

    try std.testing.expect(mat.eq(f64, A, B, 1e-12));
}

test "matrix_clone_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();

    var B = try A.clone(); defer B.deinit();

    try std.testing.expect(mat.eq(f64, A, B, 1e-12));
}

test "matrix_complex_1x1" {
    var A = try Matrix(f64         ).init(1, 1, std.testing.allocator); defer A.deinit();
    var R = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1;

    R.ptr(0, 0).* = Complex(f64).init(1, 0);

    var B = try A.complex(); defer B.deinit();

    try std.testing.expect(mat.ceq(f64, B, R, 1e-12));
}

test "matrix_complex_2x2" {
    var A = try Matrix(f64         ).init(2, 2, std.testing.allocator); defer A.deinit();
    var R = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    R.ptr(0, 0).* = Complex(f64).init(1, 0); R.ptr(0, 1).* = Complex(f64).init(2, 0);
    R.ptr(1, 0).* = Complex(f64).init(3, 0); R.ptr(1, 1).* = Complex(f64).init(4, 0);

    var B = try A.complex(); defer B.deinit();

    try std.testing.expect(mat.ceq(f64, B, R, 1e-12));
}

test "matrix_at_ptr_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();

    A.ptr(0, 0).* = 1;

    try std.testing.expect(A.at(0, 0) == 1);
}

test "matrix_at_ptr_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    try std.testing.expect(A.at(0, 0) == 1 and A.at(0, 1) == 2 and A.at(1, 0) == 3 and A.at(1, 1) == 4);
}

test "matrix_rowptr_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1;

    R.ptr(0, 0).* = 1;

    const B = A.rowptr(0);

    try std.testing.expect(mat.eq(f64, B, R, 1e-12));
}

test "matrix_rowptr_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(1, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    R.ptr(0, 0).* = 1; R.ptr(0, 1).* = 2;

    const B = A.rowptr(0);

    try std.testing.expect(mat.eq(f64, B, R, 1e-12));
}

test "matrix_vectorptr_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1;

    R.ptr(0, 0).* = 1;

    const B = A.vectorptr().matrixptr();

    try std.testing.expect(mat.eq(f64, B, R, 1e-12));
}

test "matrix_vectorptr_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(4, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    R.ptr(0, 0).* = 1; R.ptr(1, 0).* = 2; R.ptr(2, 0).* = 3; R.ptr(3, 0).* = 4;

    const B = A.vectorptr().matrixptr();

    try std.testing.expect(mat.eq(f64, B, R, 1e-12));
}

test "matrix_fill_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.fill(1);

    R.ptr(0, 0).* = 1;

    try std.testing.expect(mat.eq(f64, A, R, 1e-12));
}

test "matrix_fill_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();

    A.fill(1);

    R.ptr(0, 0).* = 1; R.ptr(0, 1).* = 1;
    R.ptr(1, 0).* = 1; R.ptr(1, 1).* = 1;

    try std.testing.expect(mat.eq(f64, A, R, 1e-12));
}

test "matrix_identity_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.identity();

    R.ptr(0, 0).* = 1;

    try std.testing.expect(mat.eq(f64, A, R, 1e-12));
}

test "matrix_identity_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();

    A.identity();

    R.ptr(0, 0).* = 1; R.ptr(0, 1).* = 0;
    R.ptr(1, 0).* = 0; R.ptr(1, 1).* = 1;

    try std.testing.expect(mat.eq(f64, A, R, 1e-12));
}

test "matrix_linspace_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.linspace(1, 1);

    R.ptr(0, 0).* = 1;

    try std.testing.expect(mat.eq(f64, A, R, 1e-12));
}

test "matrix_linspace_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();

    A.linspace(1, 4);

    R.ptr(0, 0).* = 1; R.ptr(0, 1).* = 2;
    R.ptr(1, 0).* = 3; R.ptr(1, 1).* = 4;

    try std.testing.expect(mat.eq(f64, A, R, 1e-12));
}

test "eigh_1x1" {
    var T1 = try Matrix(f64).init(1, 1, std.testing.allocator); defer T1.deinit();
    var T2 = try Matrix(f64).init(1, 1, std.testing.allocator); defer T2.deinit();

    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var J = try Matrix(f64).init(1, 1, std.testing.allocator); defer J.deinit();
    var C = try Matrix(f64).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();
    var S = try Matrix(f64).init(1, 1, std.testing.allocator); defer S.deinit();

    A.ptr(0, 0).* = 1;

    R.ptr(0, 0).* = 1;

    S.ptr(0, 0).* = 1;

    mat.eigh(f64, &J, &C, A, &T1, &T2);

    try std.testing.expect(mat.eq(f64, J, R, 1e-12) and mat.eq(f64, C, S, 1e-12));
}

test "eigh_2x2" {
    var T1 = try Matrix(f64).init(2, 2, std.testing.allocator); defer T1.deinit();
    var T2 = try Matrix(f64).init(2, 2, std.testing.allocator); defer T2.deinit();

    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var J = try Matrix(f64).init(2, 2, std.testing.allocator); defer J.deinit();
    var C = try Matrix(f64).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();
    var S = try Matrix(f64).init(2, 2, std.testing.allocator); defer S.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 2; A.ptr(1, 1).* = 1;

    R.ptr(0, 0).* = -1; R.ptr(1, 1).* = 3;

    S.ptr(0, 0).* = -std.math.sqrt1_2; S.ptr(0, 1).* = std.math.sqrt1_2;
    S.ptr(1, 0).* =  std.math.sqrt1_2; S.ptr(1, 1).* = std.math.sqrt1_2;

    mat.eigh(f64, &J, &C, A, &T1, &T2);

    try std.testing.expect(mat.eq(f64, J, R, 1e-12) and mat.eq(f64, C, S, 1e-12));
}

test "eq_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(1, 1, std.testing.allocator); defer B.deinit();

    A.ptr(0, 0).* = 1;

    B.ptr(0, 0).* = 1;

    try std.testing.expect(mat.eq(f64, A, B, 1e-12));
}

test "eq_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    B.ptr(0, 0).* = 1; B.ptr(0, 1).* = 2;
    B.ptr(1, 0).* = 3; B.ptr(1, 1).* = 4;

    try std.testing.expect(mat.eq(f64, A, B, 1e-12));
}

test "ceq_1x1" {
    var A = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer B.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2);

    B.ptr(0, 0).* = Complex(f64).init(1, 2);

    try std.testing.expect(mat.ceq(f64, A, B, 1e-12));
}

test "ceq_2x2" {
    var A = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer B.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2); A.ptr(0, 1).* = Complex(f64).init(3, 4);
    A.ptr(1, 0).* = Complex(f64).init(5, 6); A.ptr(1, 1).* = Complex(f64).init(7, 8);

    B.ptr(0, 0).* = Complex(f64).init(1, 2); B.ptr(0, 1).* = Complex(f64).init(3, 4);
    B.ptr(1, 0).* = Complex(f64).init(5, 6); B.ptr(1, 1).* = Complex(f64).init(7, 8);

    try std.testing.expect(mat.ceq(f64, A, B, 1e-12));
}

test "hjoin_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(1, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(1, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1;

    B.ptr(0, 0).* = 2;

    R.ptr(0, 0).* = 1; R.ptr(0, 1).* = 2;

    mat.hjoin(f64, &C, A, B);

    try std.testing.expect(mat.eq(f64, C, R, 1e-12));
}

test "hjoin_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(2, 4, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(2, 4, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    B.ptr(0, 0).* = 5; B.ptr(0, 1).* = 6;
    B.ptr(1, 0).* = 7; B.ptr(1, 1).* = 8;

    R.ptr(0, 0).* = 1; R.ptr(0, 1).* = 2; R.ptr(0, 2).* = 5; R.ptr(0, 3).* = 6;
    R.ptr(1, 0).* = 3; R.ptr(1, 1).* = 4; R.ptr(1, 2).* = 7; R.ptr(1, 3).* = 8;

    mat.hjoin(f64, &C, A, B);

    try std.testing.expect(mat.eq(f64, C, R, 1e-12));
}

test "mm_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 2;

    B.ptr(0, 0).* = 2;

    R.ptr(0, 0).* = 4;

    mat.mm(f64, &C, A, B);

    try std.testing.expect(mat.eq(f64, C, R, 1e-12));
}

test "mm_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    B.ptr(0, 0).* = 5; B.ptr(0, 1).* = 6;
    B.ptr(1, 0).* = 7; B.ptr(1, 1).* = 8;

    R.ptr(0, 0).* = 19; R.ptr(0, 1).* = 22;
    R.ptr(1, 0).* = 43; R.ptr(1, 1).* = 50;

    mat.mm(f64, &C, A, B);

    try std.testing.expect(mat.eq(f64, C, R, 1e-12));
}

test "cmm_1x1" {
    var A = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2);

    B.ptr(0, 0).* = Complex(f64).init(3, 4);

    R.ptr(0, 0).* = Complex(f64).init(-5, 10);

    mat.cmm(f64, &C, A, B);

    try std.testing.expect(mat.ceq(f64, C, R, 1e-12));
}

test "cmm_2x2" {
    var A = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2); A.ptr(0, 1).* = Complex(f64).init(3, 4);
    A.ptr(1, 0).* = Complex(f64).init(5, 6); A.ptr(1, 1).* = Complex(f64).init(7, 8);

    B.ptr(0, 0).* = Complex(f64).init( 9, 10); B.ptr(0, 1).* = Complex(f64).init(11, 12);
    B.ptr(1, 0).* = Complex(f64).init(13, 14); B.ptr(1, 1).* = Complex(f64).init(15, 16);

    R.ptr(0, 0).* = Complex(f64).init(-28, 122); R.ptr(0, 1).* = Complex(f64).init(-32, 142);
    R.ptr(1, 0).* = Complex(f64).init(-36, 306); R.ptr(1, 1).* = Complex(f64).init(-40, 358);

    mat.cmm(f64, &C, A, B);

    try std.testing.expect(mat.ceq(f64, C, R, 1e-12));
}

test "mam_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 2;

    B.ptr(0, 0).* = 2;

    R.ptr(0, 0).* = 4;

    mat.mam(f64, &C, A, B);

    try std.testing.expect(mat.eq(f64, C, R, 1e-12));
}

test "mam_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    B.ptr(0, 0).* = 5; B.ptr(0, 1).* = 6;
    B.ptr(1, 0).* = 7; B.ptr(1, 1).* = 8;

    R.ptr(0, 0).* = 26; R.ptr(0, 1).* = 30;
    R.ptr(1, 0).* = 38; R.ptr(1, 1).* = 44;

    mat.mam(f64, &C, A, B);

    try std.testing.expect(mat.eq(f64, C, R, 1e-12));
}

test "cmam_1x1" {
    var A = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2);

    B.ptr(0, 0).* = Complex(f64).init(3, 4);

    R.ptr(0, 0).* = Complex(f64).init(11, -2);

    mat.cmam(f64, &C, A, B);

    try std.testing.expect(mat.ceq(f64, C, R, 1e-12));
}

test "cmam_2x2" {
    var A = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2); A.ptr(0, 1).* = Complex(f64).init(3, 4);
    A.ptr(1, 0).* = Complex(f64).init(5, 6); A.ptr(1, 1).* = Complex(f64).init(7, 8);

    B.ptr(0, 0).* = Complex(f64).init( 9, 10); B.ptr(0, 1).* = Complex(f64).init(11, 12);
    B.ptr(1, 0).* = Complex(f64).init(13, 14); B.ptr(1, 1).* = Complex(f64).init(15, 16);

    R.ptr(0, 0).* = Complex(f64).init(178, -16); R.ptr(0, 1).* = Complex(f64).init(206, -20);
    R.ptr(1, 0).* = Complex(f64).init(270, -12); R.ptr(1, 1).* = Complex(f64).init(314, -16);

    mat.cmam(f64, &C, A, B);

    try std.testing.expect(mat.ceq(f64, C, R, 1e-12));
}

test "mma_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 2;

    B.ptr(0, 0).* = 2;

    R.ptr(0, 0).* = 4;

    mat.mma(f64, &C, A, B);

    try std.testing.expect(mat.eq(f64, C, R, 1e-12));
}

test "mma_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    B.ptr(0, 0).* = 5; B.ptr(0, 1).* = 6;
    B.ptr(1, 0).* = 7; B.ptr(1, 1).* = 8;

    R.ptr(0, 0).* = 17; R.ptr(0, 1).* = 23;
    R.ptr(1, 0).* = 39; R.ptr(1, 1).* = 53;

    mat.mma(f64, &C, A, B);

    try std.testing.expect(mat.eq(f64, C, R, 1e-12));
}

test "cmma_1x1" {
    var A = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2);

    B.ptr(0, 0).* = Complex(f64).init(3, 4);

    R.ptr(0, 0).* = Complex(f64).init(11, 2);

    mat.cmma(f64, &C, A, B);

    try std.testing.expect(mat.ceq(f64, C, R, 1e-12));
}

test "cmma_2x2" {
    var A = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2); A.ptr(0, 1).* = Complex(f64).init(3, 4);
    A.ptr(1, 0).* = Complex(f64).init(5, 6); A.ptr(1, 1).* = Complex(f64).init(7, 8);

    B.ptr(0, 0).* = Complex(f64).init( 9, 10); B.ptr(0, 1).* = Complex(f64).init(11, 12);
    B.ptr(1, 0).* = Complex(f64).init(13, 14); B.ptr(1, 1).* = Complex(f64).init(15, 16);

    R.ptr(0, 0).* = Complex(f64).init(110, 16); R.ptr(0, 1).* = Complex(f64).init(150, 24);
    R.ptr(1, 0).* = Complex(f64).init(278,  8); R.ptr(1, 1).* = Complex(f64).init(382, 16);

    mat.cmma(f64, &C, A, B);

    try std.testing.expect(mat.ceq(f64, C, R, 1e-12));
}

test "tmamt_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 2;

    B.ptr(0, 0).* = 2;

    R.ptr(0, 0).* = 4;

    mat.tmamt(f64, &C, A, B);

    try std.testing.expect(mat.eq(f64, C, R, 1e-12));
}

test "tmamt_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    B.ptr(0, 0).* = 5; B.ptr(0, 1).* = 6;
    B.ptr(1, 0).* = 7; B.ptr(1, 1).* = 8;

    R.ptr(0, 0).* = 23; R.ptr(0, 1).* = 34;
    R.ptr(1, 0).* = 31; R.ptr(1, 1).* = 46;

    mat.tmamt(f64, &C, A, B);

    try std.testing.expect(mat.eq(f64, C, R, 1e-12));
}

test "ctmamt_1x1" {
    var A = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2);

    B.ptr(0, 0).* = Complex(f64).init(3, 4);

    R.ptr(0, 0).* = Complex(f64).init(11, -2);

    mat.ctmamt(f64, &C, A, B);

    try std.testing.expect(mat.ceq(f64, C, R, 1e-12));
}

test "ctmamt_2x2" {
    var A = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2); A.ptr(0, 1).* = Complex(f64).init(3, 4);
    A.ptr(1, 0).* = Complex(f64).init(5, 6); A.ptr(1, 1).* = Complex(f64).init(7, 8);

    B.ptr(0, 0).* = Complex(f64).init( 9, 10); B.ptr(0, 1).* = Complex(f64).init(11, 12);
    B.ptr(1, 0).* = Complex(f64).init(13, 14); B.ptr(1, 1).* = Complex(f64).init(15, 16);

    R.ptr(0, 0).* = Complex(f64).init(156, -14); R.ptr(0, 1).* = Complex(f64).init(240, -10);
    R.ptr(1, 0).* = Complex(f64).init(212, -22); R.ptr(1, 1).* = Complex(f64).init(328, -18);

    mat.ctmamt(f64, &C, A, B);

    try std.testing.expect(mat.ceq(f64, C, R, 1e-12));
}

// FOURIER TRANSFORM ===================================================================================================================================================================================

test "fft_1d_1" {
    var A = try Vector(Complex(f64)).init(1, std.testing.allocator); defer A.deinit();
    var B = try Vector(Complex(f64)).init(1, std.testing.allocator); defer B.deinit();
    var R = try Vector(Complex(f64)).init(1, std.testing.allocator); defer R.deinit();

    A.ptr(0).* = Complex(f64).init(7, 0);

    R.ptr(0).* = Complex(f64).init(7, 0);

    ftr.fft(f64, B.data, A.data, -1);

    try std.testing.expect(vec.ceq(f64, B, R, 1e-12));
}

test "ifft_1d_1" {
    var A = try Vector(Complex(f64)).init(1, std.testing.allocator); defer A.deinit();
    var B = try Vector(Complex(f64)).init(1, std.testing.allocator); defer B.deinit();
    var R = try Vector(Complex(f64)).init(1, std.testing.allocator); defer R.deinit();

    A.ptr(0).* = Complex(f64).init(7, 0);

    R.ptr(0).* = Complex(f64).init(7, 0);

    ftr.fft(f64, B.data, A.data, 1);

    try std.testing.expect(vec.ceq(f64, B, R, 1e-12));
}

test "fft_1d_2" {
    var A = try Vector(Complex(f64)).init(2, std.testing.allocator); defer A.deinit();
    var B = try Vector(Complex(f64)).init(2, std.testing.allocator); defer B.deinit();
    var R = try Vector(Complex(f64)).init(2, std.testing.allocator); defer R.deinit();

    A.ptr(0).* = Complex(f64).init(3, 0);
    A.ptr(1).* = Complex(f64).init(1, 0);

    R.ptr(0).* = Complex(f64).init(4, 0);
    R.ptr(1).* = Complex(f64).init(2, 0);

    ftr.fft(f64, B.data, A.data, -1);

    try std.testing.expect(vec.ceq(f64, B, R, 1e-12));
}

test "ifft_1d_2" {
    var A = try Vector(Complex(f64)).init(2, std.testing.allocator); defer A.deinit();
    var B = try Vector(Complex(f64)).init(2, std.testing.allocator); defer B.deinit();
    var R = try Vector(Complex(f64)).init(2, std.testing.allocator); defer R.deinit();

    A.ptr(0).* = Complex(f64).init(4, 0);
    A.ptr(1).* = Complex(f64).init(2, 0);

    R.ptr(0).* = Complex(f64).init(3, 0);
    R.ptr(1).* = Complex(f64).init(1, 0);

    ftr.fft(f64, B.data, A.data, 1);

    try std.testing.expect(vec.ceq(f64, B, R, 1e-12));
}

test "fft_1d_4" {
    var A = try Vector(Complex(f64)).init(4, std.testing.allocator); defer A.deinit();
    var B = try Vector(Complex(f64)).init(4, std.testing.allocator); defer B.deinit();
    var R = try Vector(Complex(f64)).init(4, std.testing.allocator); defer R.deinit();

    A.ptr(0).* = Complex(f64).init(4, 0);
    A.ptr(1).* = Complex(f64).init(0, 0);
    A.ptr(2).* = Complex(f64).init(1, 0);
    A.ptr(3).* = Complex(f64).init(3, 0);

    R.ptr(0).* = Complex(f64).init(8,  0);
    R.ptr(1).* = Complex(f64).init(3,  3);
    R.ptr(2).* = Complex(f64).init(2,  0);
    R.ptr(3).* = Complex(f64).init(3, -3);

    ftr.fft(f64, B.data, A.data, -1);

    try std.testing.expect(vec.ceq(f64, B, R, 1e-12));
}

test "ifft_1d_4" {
    var A = try Vector(Complex(f64)).init(4, std.testing.allocator); defer A.deinit();
    var B = try Vector(Complex(f64)).init(4, std.testing.allocator); defer B.deinit();
    var R = try Vector(Complex(f64)).init(4, std.testing.allocator); defer R.deinit();

    A.ptr(0).* = Complex(f64).init(8,  0);
    A.ptr(1).* = Complex(f64).init(3,  3);
    A.ptr(2).* = Complex(f64).init(2,  0);
    A.ptr(3).* = Complex(f64).init(3, -3);

    R.ptr(0).* = Complex(f64).init(4, 0);
    R.ptr(1).* = Complex(f64).init(0, 0);
    R.ptr(2).* = Complex(f64).init(1, 0);
    R.ptr(3).* = Complex(f64).init(3, 0);

    ftr.fft(f64, B.data, A.data, 1);

    try std.testing.expect(vec.ceq(f64, B, R, 1e-12));
}

// QUANTUM DYNAMICS ====================================================================================================================================================================================

test "qdyn_doubleState1D_1" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.68939163804408; P.ptr(0, 1).* = 0.00018765157757;
    P.ptr(1, 0).* = 0.00018765157757; P.ptr(1, 1).* = 0.31060836195428;

    const opt = QDO(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "doubleState1D_1"
    };

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_doubleState1D_2" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.96070275892968; P.ptr(0, 1).* = 0.01414658057574;
    P.ptr(1, 0).* = 0.01414658057574; P.ptr(1, 1).* = 0.03929724106854;

    const opt = QDO(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "doubleState1D_2"
    };

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_tripleState1D_1" {
    var P = try Matrix(f64).init(3, 3, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.09406137151353; P.ptr(0, 1).* = 0.02243393031679; P.ptr(0, 2).* = 0.00177590480664;
    P.ptr(1, 0).* = 0.02243393031679; P.ptr(1, 1).* = 0.09966537676404; P.ptr(1, 2).* = 0.03403002401355;
    P.ptr(2, 0).* = 0.00177590480664; P.ptr(2, 1).* = 0.03403002401355; P.ptr(2, 2).* = 0.80627325172103;

    const opt = QDO(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 2
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "tripleState1D_1"
    };

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_tripleState1D_2" {
    var P = try Matrix(f64).init(3, 3, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.80243373396786; P.ptr(0, 1).* = 0.09151956326301; P.ptr(0, 2).* = 0.02478818099951;
    P.ptr(1, 0).* = 0.09151956326301; P.ptr(1, 1).* = 0.04418760360549; P.ptr(1, 2).* = 0.04060912220532;
    P.ptr(2, 0).* = 0.02478818099951; P.ptr(2, 1).* = 0.04060912220532; P.ptr(2, 2).* = 0.15337866242496;

    const opt = QDO(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 2
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "tripleState1D_2"
    };

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_tripleState1D_3" {
    var P = try Matrix(f64).init(3, 3, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.72007774176220; P.ptr(0, 1).* = 0.00000000137970; P.ptr(0, 2).* = 0.00000010527634;
    P.ptr(1, 0).* = 0.00000000137970; P.ptr(1, 1).* = 0.11247117674074; P.ptr(1, 2).* = 0.00000000702129;
    P.ptr(2, 0).* = 0.00000010527634; P.ptr(2, 1).* = 0.00000000702129; P.ptr(2, 2).* = 0.16745108149516;

    const opt = QDO(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "tripleState1D_3"
    };

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_tully1D_1" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.41050635034636; P.ptr(0, 1).* = 0.03954870766906;
    P.ptr(1, 0).* = 0.03954870766906; P.ptr(1, 1).* = 0.58949364965213;

    const opt = QDO(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "tully1D_1"
    };

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_tully1D_2" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.02483440861074; P.ptr(0, 1).* = 0.00092889832192;
    P.ptr(1, 0).* = 0.00092889832192; P.ptr(1, 1).* = 0.97516559138780;

    const opt = QDO(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "tully1D_2"
    };

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_tully1D_3" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.76438466406603; P.ptr(0, 1).* = 0.17171831933474;
    P.ptr(1, 0).* = 0.17171831933474; P.ptr(1, 1).* = 0.23561533593184;

    const opt = QDO(f64){
        .adiabatic = true,
        .imaginary = false,
        .iterations = 300,
        .time_step = 10,
        .grid = .{
            .limits = &[_]f64{-16, 32},
            .points = 512
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &[_]f64{15},
            .position = &[_]f64{-10},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 50
        },
        .write = .{
            .kinetic_energy = null,
            .momentum = null,
            .population = null,
            .position = null,
            .potential_energy = null,
            .total_energy = null
        },
        .potential = "tully1D_3"
    };

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

// CLASSICAL DYNAMICS ==================================================================================================================================================================================

test "cdyn_fssh_doubleState1D_1" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.72; pop.ptr(1).* = 0.28;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "doubleState1D_1",
        .type = "fssh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_doubleState1D_1" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.67; pop.ptr(1).* = 0.33;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "doubleState1D_1",
        .type = "lzsh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_fssh_doubleState1D_2" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.96; pop.ptr(1).* = 0.04;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "doubleState1D_2",
        .type = "fssh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_doubleState1D_2" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.92; pop.ptr(1).* = 0.08;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "doubleState1D_2",
        .type = "lzsh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_fssh_tripleState1D_1" {
    var pop = try Vector(f64).init(3, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.1; pop.ptr(1).* = 0.11; pop.ptr(2).* = 0.79;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 2
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "tripleState1D_1",
        .type = "fssh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_tripleState1D_1" {
    var pop = try Vector(f64).init(3, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.08; pop.ptr(1).* = 0.03; pop.ptr(2).* = 0.89;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 2
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "tripleState1D_1",
        .type = "lzsh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_fssh_tripleState1D_2" {
    var pop = try Vector(f64).init(3, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.74; pop.ptr(1).* = 0.06; pop.ptr(2).* = 0.2;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 2
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "tripleState1D_2",
        .type = "fssh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_tripleState1D_2" {
    var pop = try Vector(f64).init(3, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.79; pop.ptr(1).* = 0.13; pop.ptr(2).* = 0.08;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 2
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "tripleState1D_2",
        .type = "lzsh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_fssh_tripleState1D_3" {
    var pop = try Vector(f64).init(3, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.14; pop.ptr(1).* = 0.82; pop.ptr(2).* = 0.04;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "tripleState1D_3",
        .type = "fssh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_tripleState1D_3" {
    var pop = try Vector(f64).init(3, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.12; pop.ptr(1).* = 0.83; pop.ptr(2).* = 0.05;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "tripleState1D_3",
        .type = "lzsh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_fssh_tully1D_1" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.41; pop.ptr(1).* = 0.59;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "tully1D_1",
        .type = "fssh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_tully1D_1" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.51; pop.ptr(1).* = 0.49;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "tully1D_1",
        .type = "lzsh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_fssh_tully1D_2" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.15; pop.ptr(1).* = 0.85;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "tully1D_2",
        .type = "fssh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_tully1D_2" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.13; pop.ptr(1).* = 0.87;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "tully1D_2",
        .type = "lzsh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_fssh_tully1D_3" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0.58; pop.ptr(1).* = 0.42;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "tully1D_3",
        .type = "fssh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_tully1D_3" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    pop.ptr(0).* = 0; pop.ptr(1).* = 1;

    const opt = CDO(f64){
        .adiabatic = true,
        .derivative_step = 0.001,
        .iterations = 3000,
        .seed = 1,
        .time_step = 1,
        .trajectories = 100,
        .initial_conditions = .{
            .mass = 2000,
            .momentum_mean = &[_]f64{15},
            .momentum_std = &[_]f64{1},
            .position_mean = &[_]f64{-10},
            .position_std = &[_]f64{0.5},
            .state = 1
        },
        .log_intervals = .{
            .iteration = 500,
            .trajectory = 100
        },
        .write = .{
            .fssh_coefficient_mean = null,
            .kinetic_energy_mean = null,
            .momentum_mean = null,
            .population_mean = null,
            .position_mean = null,
            .potential_energy_mean = null,
            .total_energy_mean = null
        },
        .potential = "tully1D_3",
        .type = "lzsh"
    };

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}
