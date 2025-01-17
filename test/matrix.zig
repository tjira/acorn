const std = @import("std");

const mat = @import("acorn").mat;

const Matrix = @import("acorn").Matrix;

const allocator = std.testing.allocator;

test "mat_init" {
    const A = try mat.Matrix(f64).init(1   , 1   , allocator); defer A.deinit();
    const B = try mat.Matrix(f64).init(1   , 10  , allocator); defer B.deinit();
    const C = try mat.Matrix(f64).init(1   , 100 , allocator); defer C.deinit();
    const D = try mat.Matrix(f64).init(1   , 1000, allocator); defer D.deinit();
    const E = try mat.Matrix(f64).init(10  , 1   , allocator); defer E.deinit();
    const F = try mat.Matrix(f64).init(10  , 10  , allocator); defer F.deinit();
    const G = try mat.Matrix(f64).init(10  , 100 , allocator); defer G.deinit();
    const H = try mat.Matrix(f64).init(10  , 1000, allocator); defer H.deinit();
    const I = try mat.Matrix(f64).init(100 , 1   , allocator); defer I.deinit();
    const J = try mat.Matrix(f64).init(100 , 10  , allocator); defer J.deinit();
    const K = try mat.Matrix(f64).init(100 , 100 , allocator); defer K.deinit();
    const L = try mat.Matrix(f64).init(100 , 1000, allocator); defer L.deinit();
    const M = try mat.Matrix(f64).init(1000, 1   , allocator); defer M.deinit();
    const N = try mat.Matrix(f64).init(1000, 10  , allocator); defer N.deinit();
    const O = try mat.Matrix(f64).init(1000, 100 , allocator); defer O.deinit();
    const P = try mat.Matrix(f64).init(1000, 1000, allocator); defer P.deinit();

    try std.testing.expect(A.rows == 1    and A.cols == 1   );
    try std.testing.expect(B.rows == 1    and B.cols == 10  );
    try std.testing.expect(C.rows == 1    and C.cols == 100 );
    try std.testing.expect(D.rows == 1    and D.cols == 1000);
    try std.testing.expect(E.rows == 10   and E.cols == 1   );
    try std.testing.expect(F.rows == 10   and F.cols == 10  );
    try std.testing.expect(G.rows == 10   and G.cols == 100 );
    try std.testing.expect(H.rows == 10   and H.cols == 1000);
    try std.testing.expect(I.rows == 100  and I.cols == 1   );
    try std.testing.expect(J.rows == 100  and J.cols == 10  );
    try std.testing.expect(K.rows == 100  and K.cols == 100 );
    try std.testing.expect(L.rows == 100  and L.cols == 1000);
    try std.testing.expect(M.rows == 1000 and M.cols == 1   );
    try std.testing.expect(N.rows == 1000 and N.cols == 10  );
    try std.testing.expect(O.rows == 1000 and O.cols == 100 );
    try std.testing.expect(P.rows == 1000 and P.cols == 1000);
}
