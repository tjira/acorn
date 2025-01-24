const std = @import("std");

const helper = @import("acorn").helper;
const matrix = @import("acorn").matrix;

const Matrix = @import("acorn").Matrix;

const allocator = std.testing.allocator;

test "matrix_at" {
    const A = try matrix.Matrix(f64).init(1   , 1   , allocator); defer A.deinit();
    const B = try matrix.Matrix(f64).init(1   , 10  , allocator); defer B.deinit();
    const C = try matrix.Matrix(f64).init(1   , 100 , allocator); defer C.deinit();
    const D = try matrix.Matrix(f64).init(1   , 1000, allocator); defer D.deinit();
    const E = try matrix.Matrix(f64).init(10  , 1   , allocator); defer E.deinit();
    const F = try matrix.Matrix(f64).init(10  , 10  , allocator); defer F.deinit();
    const G = try matrix.Matrix(f64).init(10  , 100 , allocator); defer G.deinit();
    const H = try matrix.Matrix(f64).init(10  , 1000, allocator); defer H.deinit();
    const I = try matrix.Matrix(f64).init(100 , 1   , allocator); defer I.deinit();
    const J = try matrix.Matrix(f64).init(100 , 10  , allocator); defer J.deinit();
    const K = try matrix.Matrix(f64).init(100 , 100 , allocator); defer K.deinit();
    const L = try matrix.Matrix(f64).init(100 , 1000, allocator); defer L.deinit();
    const M = try matrix.Matrix(f64).init(1000, 1   , allocator); defer M.deinit();
    const N = try matrix.Matrix(f64).init(1000, 10  , allocator); defer N.deinit();
    const O = try matrix.Matrix(f64).init(1000, 100 , allocator); defer O.deinit();
    const P = try matrix.Matrix(f64).init(1000, 1000, allocator); defer P.deinit();

    for (0..A.rows) |i| for (0..A.cols) |j| {A.ptr(i, j).* = helper.asfloat(f64, i * A.cols + j);};
    for (0..B.rows) |i| for (0..B.cols) |j| {B.ptr(i, j).* = helper.asfloat(f64, i * B.cols + j);};
    for (0..C.rows) |i| for (0..C.cols) |j| {C.ptr(i, j).* = helper.asfloat(f64, i * C.cols + j);};
    for (0..D.rows) |i| for (0..D.cols) |j| {D.ptr(i, j).* = helper.asfloat(f64, i * D.cols + j);};
    for (0..E.rows) |i| for (0..E.cols) |j| {E.ptr(i, j).* = helper.asfloat(f64, i * E.cols + j);};
    for (0..F.rows) |i| for (0..F.cols) |j| {F.ptr(i, j).* = helper.asfloat(f64, i * F.cols + j);};
    for (0..G.rows) |i| for (0..G.cols) |j| {G.ptr(i, j).* = helper.asfloat(f64, i * G.cols + j);};
    for (0..H.rows) |i| for (0..H.cols) |j| {H.ptr(i, j).* = helper.asfloat(f64, i * H.cols + j);};
    for (0..I.rows) |i| for (0..I.cols) |j| {I.ptr(i, j).* = helper.asfloat(f64, i * I.cols + j);};
    for (0..J.rows) |i| for (0..J.cols) |j| {J.ptr(i, j).* = helper.asfloat(f64, i * J.cols + j);};
    for (0..K.rows) |i| for (0..K.cols) |j| {K.ptr(i, j).* = helper.asfloat(f64, i * K.cols + j);};
    for (0..L.rows) |i| for (0..L.cols) |j| {L.ptr(i, j).* = helper.asfloat(f64, i * L.cols + j);};
    for (0..M.rows) |i| for (0..M.cols) |j| {M.ptr(i, j).* = helper.asfloat(f64, i * M.cols + j);};
    for (0..N.rows) |i| for (0..N.cols) |j| {N.ptr(i, j).* = helper.asfloat(f64, i * N.cols + j);};
    for (0..O.rows) |i| for (0..O.cols) |j| {O.ptr(i, j).* = helper.asfloat(f64, i * O.cols + j);};
    for (0..P.rows) |i| for (0..P.cols) |j| {P.ptr(i, j).* = helper.asfloat(f64, i * P.cols + j);};

    for (0..A.rows) |i| for (0..A.cols) |j| {try std.testing.expect(A.at(i, j) == helper.asfloat(f64, i * A.cols + j));};
    for (0..B.rows) |i| for (0..B.cols) |j| {try std.testing.expect(B.at(i, j) == helper.asfloat(f64, i * B.cols + j));};
    for (0..C.rows) |i| for (0..C.cols) |j| {try std.testing.expect(C.at(i, j) == helper.asfloat(f64, i * C.cols + j));};
    for (0..D.rows) |i| for (0..D.cols) |j| {try std.testing.expect(D.at(i, j) == helper.asfloat(f64, i * D.cols + j));};
    for (0..E.rows) |i| for (0..E.cols) |j| {try std.testing.expect(E.at(i, j) == helper.asfloat(f64, i * E.cols + j));};
    for (0..F.rows) |i| for (0..F.cols) |j| {try std.testing.expect(F.at(i, j) == helper.asfloat(f64, i * F.cols + j));};
    for (0..G.rows) |i| for (0..G.cols) |j| {try std.testing.expect(G.at(i, j) == helper.asfloat(f64, i * G.cols + j));};
    for (0..H.rows) |i| for (0..H.cols) |j| {try std.testing.expect(H.at(i, j) == helper.asfloat(f64, i * H.cols + j));};
    for (0..I.rows) |i| for (0..I.cols) |j| {try std.testing.expect(I.at(i, j) == helper.asfloat(f64, i * I.cols + j));};
    for (0..J.rows) |i| for (0..J.cols) |j| {try std.testing.expect(J.at(i, j) == helper.asfloat(f64, i * J.cols + j));};
    for (0..K.rows) |i| for (0..K.cols) |j| {try std.testing.expect(K.at(i, j) == helper.asfloat(f64, i * K.cols + j));};
    for (0..L.rows) |i| for (0..L.cols) |j| {try std.testing.expect(L.at(i, j) == helper.asfloat(f64, i * L.cols + j));};
    for (0..M.rows) |i| for (0..M.cols) |j| {try std.testing.expect(M.at(i, j) == helper.asfloat(f64, i * M.cols + j));};
    for (0..N.rows) |i| for (0..N.cols) |j| {try std.testing.expect(N.at(i, j) == helper.asfloat(f64, i * N.cols + j));};
    for (0..O.rows) |i| for (0..O.cols) |j| {try std.testing.expect(O.at(i, j) == helper.asfloat(f64, i * O.cols + j));};
    for (0..P.rows) |i| for (0..P.cols) |j| {try std.testing.expect(P.at(i, j) == helper.asfloat(f64, i * P.cols + j));};
}

test "matrix_clone" {
    const A = try matrix.Matrix(f64).init(1   , 1   , allocator); defer A.deinit();
    const B = try matrix.Matrix(f64).init(1   , 10  , allocator); defer B.deinit();
    const C = try matrix.Matrix(f64).init(1   , 100 , allocator); defer C.deinit();
    const D = try matrix.Matrix(f64).init(1   , 1000, allocator); defer D.deinit();
    const E = try matrix.Matrix(f64).init(10  , 1   , allocator); defer E.deinit();
    const F = try matrix.Matrix(f64).init(10  , 10  , allocator); defer F.deinit();
    const G = try matrix.Matrix(f64).init(10  , 100 , allocator); defer G.deinit();
    const H = try matrix.Matrix(f64).init(10  , 1000, allocator); defer H.deinit();
    const I = try matrix.Matrix(f64).init(100 , 1   , allocator); defer I.deinit();
    const J = try matrix.Matrix(f64).init(100 , 10  , allocator); defer J.deinit();
    const K = try matrix.Matrix(f64).init(100 , 100 , allocator); defer K.deinit();
    const L = try matrix.Matrix(f64).init(100 , 1000, allocator); defer L.deinit();
    const M = try matrix.Matrix(f64).init(1000, 1   , allocator); defer M.deinit();
    const N = try matrix.Matrix(f64).init(1000, 10  , allocator); defer N.deinit();
    const O = try matrix.Matrix(f64).init(1000, 100 , allocator); defer O.deinit();
    const P = try matrix.Matrix(f64).init(1000, 1000, allocator); defer P.deinit();

    A.randn(0, 1, 0); B.randn(0, 1, 1); C.randn(0, 1, 2 ); D.randn(0, 1, 3 ); E.randn(0, 1, 4 ); F.randn(0, 1, 5 ); G.randn(0, 1, 6 ); H.randn(0, 1, 7 );
    I.randn(0, 1, 8); J.randn(0, 1, 9); K.randn(0, 1, 10); L.randn(0, 1, 11); M.randn(0, 1, 12); N.randn(0, 1, 13); O.randn(0, 1, 14); P.randn(0, 1, 15);

    const AA = try A.clone(); defer AA.deinit();
    const BB = try B.clone(); defer BB.deinit();
    const CC = try C.clone(); defer CC.deinit();
    const DD = try D.clone(); defer DD.deinit();
    const EE = try E.clone(); defer EE.deinit();
    const FF = try F.clone(); defer FF.deinit();
    const GG = try G.clone(); defer GG.deinit();
    const HH = try H.clone(); defer HH.deinit();
    const II = try I.clone(); defer II.deinit();
    const JJ = try J.clone(); defer JJ.deinit();
    const KK = try K.clone(); defer KK.deinit();
    const LL = try L.clone(); defer LL.deinit();
    const MM = try M.clone(); defer MM.deinit();
    const NN = try N.clone(); defer NN.deinit();
    const OO = try O.clone(); defer OO.deinit();
    const PP = try P.clone(); defer PP.deinit();

    for (0..A.rows) |i| for (0..A.cols) |j| {try std.testing.expect(A.at(i, j) == AA.at(i, j));};
    for (0..B.rows) |i| for (0..B.cols) |j| {try std.testing.expect(B.at(i, j) == BB.at(i, j));};
    for (0..C.rows) |i| for (0..C.cols) |j| {try std.testing.expect(C.at(i, j) == CC.at(i, j));};
    for (0..D.rows) |i| for (0..D.cols) |j| {try std.testing.expect(D.at(i, j) == DD.at(i, j));};
    for (0..E.rows) |i| for (0..E.cols) |j| {try std.testing.expect(E.at(i, j) == EE.at(i, j));};
    for (0..F.rows) |i| for (0..F.cols) |j| {try std.testing.expect(F.at(i, j) == FF.at(i, j));};
    for (0..G.rows) |i| for (0..G.cols) |j| {try std.testing.expect(G.at(i, j) == GG.at(i, j));};
    for (0..H.rows) |i| for (0..H.cols) |j| {try std.testing.expect(H.at(i, j) == HH.at(i, j));};
    for (0..I.rows) |i| for (0..I.cols) |j| {try std.testing.expect(I.at(i, j) == II.at(i, j));};
    for (0..J.rows) |i| for (0..J.cols) |j| {try std.testing.expect(J.at(i, j) == JJ.at(i, j));};
    for (0..K.rows) |i| for (0..K.cols) |j| {try std.testing.expect(K.at(i, j) == KK.at(i, j));};
    for (0..L.rows) |i| for (0..L.cols) |j| {try std.testing.expect(L.at(i, j) == LL.at(i, j));};
    for (0..M.rows) |i| for (0..M.cols) |j| {try std.testing.expect(M.at(i, j) == MM.at(i, j));};
    for (0..N.rows) |i| for (0..N.cols) |j| {try std.testing.expect(N.at(i, j) == NN.at(i, j));};
    for (0..O.rows) |i| for (0..O.cols) |j| {try std.testing.expect(O.at(i, j) == OO.at(i, j));};
    for (0..P.rows) |i| for (0..P.cols) |j| {try std.testing.expect(P.at(i, j) == PP.at(i, j));};
}

test "matrix_complex" {
    const A = try matrix.Matrix(f64).init(1   , 1   , allocator); defer A.deinit();
    const B = try matrix.Matrix(f64).init(1   , 10  , allocator); defer B.deinit();
    const C = try matrix.Matrix(f64).init(1   , 100 , allocator); defer C.deinit();
    const D = try matrix.Matrix(f64).init(1   , 1000, allocator); defer D.deinit();
    const E = try matrix.Matrix(f64).init(10  , 1   , allocator); defer E.deinit();
    const F = try matrix.Matrix(f64).init(10  , 10  , allocator); defer F.deinit();
    const G = try matrix.Matrix(f64).init(10  , 100 , allocator); defer G.deinit();
    const H = try matrix.Matrix(f64).init(10  , 1000, allocator); defer H.deinit();
    const I = try matrix.Matrix(f64).init(100 , 1   , allocator); defer I.deinit();
    const J = try matrix.Matrix(f64).init(100 , 10  , allocator); defer J.deinit();
    const K = try matrix.Matrix(f64).init(100 , 100 , allocator); defer K.deinit();
    const L = try matrix.Matrix(f64).init(100 , 1000, allocator); defer L.deinit();
    const M = try matrix.Matrix(f64).init(1000, 1   , allocator); defer M.deinit();
    const N = try matrix.Matrix(f64).init(1000, 10  , allocator); defer N.deinit();
    const O = try matrix.Matrix(f64).init(1000, 100 , allocator); defer O.deinit();
    const P = try matrix.Matrix(f64).init(1000, 1000, allocator); defer P.deinit();

    A.randn(0, 1, 0); B.randn(0, 1, 1); C.randn(0, 1, 2 ); D.randn(0, 1, 3 ); E.randn(0, 1, 4 ); F.randn(0, 1, 5 ); G.randn(0, 1, 6 ); H.randn(0, 1, 7 );
    I.randn(0, 1, 8); J.randn(0, 1, 9); K.randn(0, 1, 10); L.randn(0, 1, 11); M.randn(0, 1, 12); N.randn(0, 1, 13); O.randn(0, 1, 14); P.randn(0, 1, 15);

    const AA = try A.complex(); defer AA.deinit();
    const BB = try B.complex(); defer BB.deinit();
    const CC = try C.complex(); defer CC.deinit();
    const DD = try D.complex(); defer DD.deinit();
    const EE = try E.complex(); defer EE.deinit();
    const FF = try F.complex(); defer FF.deinit();
    const GG = try G.complex(); defer GG.deinit();
    const HH = try H.complex(); defer HH.deinit();
    const II = try I.complex(); defer II.deinit();
    const JJ = try J.complex(); defer JJ.deinit();
    const KK = try K.complex(); defer KK.deinit();
    const LL = try L.complex(); defer LL.deinit();
    const MM = try M.complex(); defer MM.deinit();
    const NN = try N.complex(); defer NN.deinit();
    const OO = try O.complex(); defer OO.deinit();
    const PP = try P.complex(); defer PP.deinit();

    for (0..A.rows) |i| for (0..A.cols) |j| {try std.testing.expect(A.at(i, j) == AA.at(i, j).re and AA.at(i, j).im == 0);};
    for (0..B.rows) |i| for (0..B.cols) |j| {try std.testing.expect(B.at(i, j) == BB.at(i, j).re and BB.at(i, j).im == 0);};
    for (0..C.rows) |i| for (0..C.cols) |j| {try std.testing.expect(C.at(i, j) == CC.at(i, j).re and CC.at(i, j).im == 0);};
    for (0..D.rows) |i| for (0..D.cols) |j| {try std.testing.expect(D.at(i, j) == DD.at(i, j).re and DD.at(i, j).im == 0);};
    for (0..E.rows) |i| for (0..E.cols) |j| {try std.testing.expect(E.at(i, j) == EE.at(i, j).re and EE.at(i, j).im == 0);};
    for (0..F.rows) |i| for (0..F.cols) |j| {try std.testing.expect(F.at(i, j) == FF.at(i, j).re and FF.at(i, j).im == 0);};
    for (0..G.rows) |i| for (0..G.cols) |j| {try std.testing.expect(G.at(i, j) == GG.at(i, j).re and GG.at(i, j).im == 0);};
    for (0..H.rows) |i| for (0..H.cols) |j| {try std.testing.expect(H.at(i, j) == HH.at(i, j).re and HH.at(i, j).im == 0);};
    for (0..I.rows) |i| for (0..I.cols) |j| {try std.testing.expect(I.at(i, j) == II.at(i, j).re and II.at(i, j).im == 0);};
    for (0..J.rows) |i| for (0..J.cols) |j| {try std.testing.expect(J.at(i, j) == JJ.at(i, j).re and JJ.at(i, j).im == 0);};
    for (0..K.rows) |i| for (0..K.cols) |j| {try std.testing.expect(K.at(i, j) == KK.at(i, j).re and KK.at(i, j).im == 0);};
    for (0..L.rows) |i| for (0..L.cols) |j| {try std.testing.expect(L.at(i, j) == LL.at(i, j).re and LL.at(i, j).im == 0);};
    for (0..M.rows) |i| for (0..M.cols) |j| {try std.testing.expect(M.at(i, j) == MM.at(i, j).re and MM.at(i, j).im == 0);};
    for (0..N.rows) |i| for (0..N.cols) |j| {try std.testing.expect(N.at(i, j) == NN.at(i, j).re and NN.at(i, j).im == 0);};
    for (0..O.rows) |i| for (0..O.cols) |j| {try std.testing.expect(O.at(i, j) == OO.at(i, j).re and OO.at(i, j).im == 0);};
    for (0..P.rows) |i| for (0..P.cols) |j| {try std.testing.expect(P.at(i, j) == PP.at(i, j).re and PP.at(i, j).im == 0);};
}

test "matrix_fill" {
    const A = try matrix.Matrix(f64).init(1   , 1   , allocator); defer A.deinit();
    const B = try matrix.Matrix(f64).init(1   , 10  , allocator); defer B.deinit();
    const C = try matrix.Matrix(f64).init(1   , 100 , allocator); defer C.deinit();
    const D = try matrix.Matrix(f64).init(1   , 1000, allocator); defer D.deinit();
    const E = try matrix.Matrix(f64).init(10  , 1   , allocator); defer E.deinit();
    const F = try matrix.Matrix(f64).init(10  , 10  , allocator); defer F.deinit();
    const G = try matrix.Matrix(f64).init(10  , 100 , allocator); defer G.deinit();
    const H = try matrix.Matrix(f64).init(10  , 1000, allocator); defer H.deinit();
    const I = try matrix.Matrix(f64).init(100 , 1   , allocator); defer I.deinit();
    const J = try matrix.Matrix(f64).init(100 , 10  , allocator); defer J.deinit();
    const K = try matrix.Matrix(f64).init(100 , 100 , allocator); defer K.deinit();
    const L = try matrix.Matrix(f64).init(100 , 1000, allocator); defer L.deinit();
    const M = try matrix.Matrix(f64).init(1000, 1   , allocator); defer M.deinit();
    const N = try matrix.Matrix(f64).init(1000, 10  , allocator); defer N.deinit();
    const O = try matrix.Matrix(f64).init(1000, 100 , allocator); defer O.deinit();
    const P = try matrix.Matrix(f64).init(1000, 1000, allocator); defer P.deinit();

    A.fill(0); B.fill(1); C.fill(2 ); D.fill(3 ); E.fill(4 ); F.fill(5 ); G.fill(6 ); H.fill(7 );
    I.fill(8); J.fill(9); K.fill(10); L.fill(11); M.fill(12); N.fill(13); O.fill(14); P.fill(15);

    for (0..A.rows) |i| for (0..A.cols) |j| {try std.testing.expect(A.at(i, j) == 0 );};
    for (0..B.rows) |i| for (0..B.cols) |j| {try std.testing.expect(B.at(i, j) == 1 );};
    for (0..C.rows) |i| for (0..C.cols) |j| {try std.testing.expect(C.at(i, j) == 2 );};
    for (0..D.rows) |i| for (0..D.cols) |j| {try std.testing.expect(D.at(i, j) == 3 );};
    for (0..E.rows) |i| for (0..E.cols) |j| {try std.testing.expect(E.at(i, j) == 4 );};
    for (0..F.rows) |i| for (0..F.cols) |j| {try std.testing.expect(F.at(i, j) == 5 );};
    for (0..G.rows) |i| for (0..G.cols) |j| {try std.testing.expect(G.at(i, j) == 6 );};
    for (0..H.rows) |i| for (0..H.cols) |j| {try std.testing.expect(H.at(i, j) == 7 );};
    for (0..I.rows) |i| for (0..I.cols) |j| {try std.testing.expect(I.at(i, j) == 8 );};
    for (0..J.rows) |i| for (0..J.cols) |j| {try std.testing.expect(J.at(i, j) == 9 );};
    for (0..K.rows) |i| for (0..K.cols) |j| {try std.testing.expect(K.at(i, j) == 10);};
    for (0..L.rows) |i| for (0..L.cols) |j| {try std.testing.expect(L.at(i, j) == 11);};
    for (0..M.rows) |i| for (0..M.cols) |j| {try std.testing.expect(M.at(i, j) == 12);};
    for (0..N.rows) |i| for (0..N.cols) |j| {try std.testing.expect(N.at(i, j) == 13);};
    for (0..O.rows) |i| for (0..O.cols) |j| {try std.testing.expect(O.at(i, j) == 14);};
    for (0..P.rows) |i| for (0..P.cols) |j| {try std.testing.expect(P.at(i, j) == 15);};
}

test "matrix_identity" {
    const A = try matrix.Matrix(f64).init(1   , 1   , allocator); defer A.deinit();
    const B = try matrix.Matrix(f64).init(1   , 10  , allocator); defer B.deinit();
    const C = try matrix.Matrix(f64).init(1   , 100 , allocator); defer C.deinit();
    const D = try matrix.Matrix(f64).init(1   , 1000, allocator); defer D.deinit();
    const E = try matrix.Matrix(f64).init(10  , 1   , allocator); defer E.deinit();
    const F = try matrix.Matrix(f64).init(10  , 10  , allocator); defer F.deinit();
    const G = try matrix.Matrix(f64).init(10  , 100 , allocator); defer G.deinit();
    const H = try matrix.Matrix(f64).init(10  , 1000, allocator); defer H.deinit();
    const I = try matrix.Matrix(f64).init(100 , 1   , allocator); defer I.deinit();
    const J = try matrix.Matrix(f64).init(100 , 10  , allocator); defer J.deinit();
    const K = try matrix.Matrix(f64).init(100 , 100 , allocator); defer K.deinit();
    const L = try matrix.Matrix(f64).init(100 , 1000, allocator); defer L.deinit();
    const M = try matrix.Matrix(f64).init(1000, 1   , allocator); defer M.deinit();
    const N = try matrix.Matrix(f64).init(1000, 10  , allocator); defer N.deinit();
    const O = try matrix.Matrix(f64).init(1000, 100 , allocator); defer O.deinit();
    const P = try matrix.Matrix(f64).init(1000, 1000, allocator); defer P.deinit();

    A.identity(); B.identity(); C.identity(); D.identity(); E.identity(); F.identity(); G.identity(); H.identity();
    I.identity(); J.identity(); K.identity(); L.identity(); M.identity(); N.identity(); O.identity(); P.identity();

    for (0..A.rows) |i| for (0..A.cols) |j| {if (i == j) try std.testing.expect(A.at(i, j) == 1) else try std.testing.expect(A.at(i, j) == 0);};
    for (0..B.rows) |i| for (0..B.cols) |j| {if (i == j) try std.testing.expect(B.at(i, j) == 1) else try std.testing.expect(B.at(i, j) == 0);};
    for (0..C.rows) |i| for (0..C.cols) |j| {if (i == j) try std.testing.expect(C.at(i, j) == 1) else try std.testing.expect(C.at(i, j) == 0);};
    for (0..D.rows) |i| for (0..D.cols) |j| {if (i == j) try std.testing.expect(D.at(i, j) == 1) else try std.testing.expect(D.at(i, j) == 0);};
    for (0..E.rows) |i| for (0..E.cols) |j| {if (i == j) try std.testing.expect(E.at(i, j) == 1) else try std.testing.expect(E.at(i, j) == 0);};
    for (0..F.rows) |i| for (0..F.cols) |j| {if (i == j) try std.testing.expect(F.at(i, j) == 1) else try std.testing.expect(F.at(i, j) == 0);};
    for (0..G.rows) |i| for (0..G.cols) |j| {if (i == j) try std.testing.expect(G.at(i, j) == 1) else try std.testing.expect(G.at(i, j) == 0);};
    for (0..H.rows) |i| for (0..H.cols) |j| {if (i == j) try std.testing.expect(H.at(i, j) == 1) else try std.testing.expect(H.at(i, j) == 0);};
    for (0..I.rows) |i| for (0..I.cols) |j| {if (i == j) try std.testing.expect(I.at(i, j) == 1) else try std.testing.expect(I.at(i, j) == 0);};
    for (0..J.rows) |i| for (0..J.cols) |j| {if (i == j) try std.testing.expect(J.at(i, j) == 1) else try std.testing.expect(J.at(i, j) == 0);};
    for (0..K.rows) |i| for (0..K.cols) |j| {if (i == j) try std.testing.expect(K.at(i, j) == 1) else try std.testing.expect(K.at(i, j) == 0);};
    for (0..L.rows) |i| for (0..L.cols) |j| {if (i == j) try std.testing.expect(L.at(i, j) == 1) else try std.testing.expect(L.at(i, j) == 0);};
    for (0..M.rows) |i| for (0..M.cols) |j| {if (i == j) try std.testing.expect(M.at(i, j) == 1) else try std.testing.expect(M.at(i, j) == 0);};
    for (0..N.rows) |i| for (0..N.cols) |j| {if (i == j) try std.testing.expect(N.at(i, j) == 1) else try std.testing.expect(N.at(i, j) == 0);};
    for (0..O.rows) |i| for (0..O.cols) |j| {if (i == j) try std.testing.expect(O.at(i, j) == 1) else try std.testing.expect(O.at(i, j) == 0);};
    for (0..P.rows) |i| for (0..P.cols) |j| {if (i == j) try std.testing.expect(P.at(i, j) == 1) else try std.testing.expect(P.at(i, j) == 0);};
}

test "matrix_init" {
    const A = try matrix.Matrix(f64).init(1   , 1   , allocator); defer A.deinit();
    const B = try matrix.Matrix(f64).init(1   , 10  , allocator); defer B.deinit();
    const C = try matrix.Matrix(f64).init(1   , 100 , allocator); defer C.deinit();
    const D = try matrix.Matrix(f64).init(1   , 1000, allocator); defer D.deinit();
    const E = try matrix.Matrix(f64).init(10  , 1   , allocator); defer E.deinit();
    const F = try matrix.Matrix(f64).init(10  , 10  , allocator); defer F.deinit();
    const G = try matrix.Matrix(f64).init(10  , 100 , allocator); defer G.deinit();
    const H = try matrix.Matrix(f64).init(10  , 1000, allocator); defer H.deinit();
    const I = try matrix.Matrix(f64).init(100 , 1   , allocator); defer I.deinit();
    const J = try matrix.Matrix(f64).init(100 , 10  , allocator); defer J.deinit();
    const K = try matrix.Matrix(f64).init(100 , 100 , allocator); defer K.deinit();
    const L = try matrix.Matrix(f64).init(100 , 1000, allocator); defer L.deinit();
    const M = try matrix.Matrix(f64).init(1000, 1   , allocator); defer M.deinit();
    const N = try matrix.Matrix(f64).init(1000, 10  , allocator); defer N.deinit();
    const O = try matrix.Matrix(f64).init(1000, 100 , allocator); defer O.deinit();
    const P = try matrix.Matrix(f64).init(1000, 1000, allocator); defer P.deinit();

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

test "matrix_io" {
    const A = try matrix.Matrix(f64).init(1   , 1   , allocator); defer A.deinit();
    const B = try matrix.Matrix(f64).init(1   , 10  , allocator); defer B.deinit();
    const C = try matrix.Matrix(f64).init(1   , 100 , allocator); defer C.deinit();
    const D = try matrix.Matrix(f64).init(1   , 1000, allocator); defer D.deinit();
    const E = try matrix.Matrix(f64).init(10  , 1   , allocator); defer E.deinit();
    const F = try matrix.Matrix(f64).init(10  , 10  , allocator); defer F.deinit();
    const G = try matrix.Matrix(f64).init(10  , 100 , allocator); defer G.deinit();
    const H = try matrix.Matrix(f64).init(10  , 1000, allocator); defer H.deinit();
    const I = try matrix.Matrix(f64).init(100 , 1   , allocator); defer I.deinit();
    const J = try matrix.Matrix(f64).init(100 , 10  , allocator); defer J.deinit();
    const K = try matrix.Matrix(f64).init(100 , 100 , allocator); defer K.deinit();
    const L = try matrix.Matrix(f64).init(100 , 1000, allocator); defer L.deinit();
    const M = try matrix.Matrix(f64).init(1000, 1   , allocator); defer M.deinit();
    const N = try matrix.Matrix(f64).init(1000, 10  , allocator); defer N.deinit();
    const O = try matrix.Matrix(f64).init(1000, 100 , allocator); defer O.deinit();
    const P = try matrix.Matrix(f64).init(1000, 1000, allocator); defer P.deinit();

    A.randn(0, 1, 0); B.randn(0, 1, 1); C.randn(0, 1, 2 ); D.randn(0, 1, 3 ); E.randn(0, 1, 4 ); F.randn(0, 1, 5 ); G.randn(0, 1, 6 ); H.randn(0, 1, 7 );
    I.randn(0, 1, 8); J.randn(0, 1, 9); K.randn(0, 1, 10); L.randn(0, 1, 11); M.randn(0, 1, 12); N.randn(0, 1, 13); O.randn(0, 1, 14); P.randn(0, 1, 15);

    try A.write("A.mat"); try B.write("B.mat"); try C.write("C.mat"); try D.write("D.mat"); try E.write("E.mat"); try F.write("F.mat"); try G.write("G.mat"); try H.write("H.mat");
    try I.write("I.mat"); try J.write("J.mat"); try K.write("K.mat"); try L.write("L.mat"); try M.write("M.mat"); try N.write("N.mat"); try O.write("O.mat"); try P.write("P.mat");

    const AA = try matrix.read(f64, "A.mat", allocator); defer AA.deinit();
    const BB = try matrix.read(f64, "B.mat", allocator); defer BB.deinit();
    const CC = try matrix.read(f64, "C.mat", allocator); defer CC.deinit();
    const DD = try matrix.read(f64, "D.mat", allocator); defer DD.deinit();
    const EE = try matrix.read(f64, "E.mat", allocator); defer EE.deinit();
    const FF = try matrix.read(f64, "F.mat", allocator); defer FF.deinit();
    const GG = try matrix.read(f64, "G.mat", allocator); defer GG.deinit();
    const HH = try matrix.read(f64, "H.mat", allocator); defer HH.deinit();
    const II = try matrix.read(f64, "I.mat", allocator); defer II.deinit();
    const JJ = try matrix.read(f64, "J.mat", allocator); defer JJ.deinit();
    const KK = try matrix.read(f64, "K.mat", allocator); defer KK.deinit();
    const LL = try matrix.read(f64, "L.mat", allocator); defer LL.deinit();
    const MM = try matrix.read(f64, "M.mat", allocator); defer MM.deinit();
    const NN = try matrix.read(f64, "N.mat", allocator); defer NN.deinit();
    const OO = try matrix.read(f64, "O.mat", allocator); defer OO.deinit();
    const PP = try matrix.read(f64, "P.mat", allocator); defer PP.deinit();

    try std.testing.expect(A.rows == AA.rows and A.cols == AA.cols);
    try std.testing.expect(B.rows == BB.rows and B.cols == BB.cols);
    try std.testing.expect(C.rows == CC.rows and C.cols == CC.cols);
    try std.testing.expect(D.rows == DD.rows and D.cols == DD.cols);
    try std.testing.expect(E.rows == EE.rows and E.cols == EE.cols);
    try std.testing.expect(F.rows == FF.rows and F.cols == FF.cols);
    try std.testing.expect(G.rows == GG.rows and G.cols == GG.cols);
    try std.testing.expect(H.rows == HH.rows and H.cols == HH.cols);
    try std.testing.expect(I.rows == II.rows and I.cols == II.cols);
    try std.testing.expect(J.rows == JJ.rows and J.cols == JJ.cols);
    try std.testing.expect(K.rows == KK.rows and K.cols == KK.cols);
    try std.testing.expect(L.rows == LL.rows and L.cols == LL.cols);
    try std.testing.expect(M.rows == MM.rows and M.cols == MM.cols);
    try std.testing.expect(N.rows == NN.rows and N.cols == NN.cols);
    try std.testing.expect(O.rows == OO.rows and O.cols == OO.cols);
    try std.testing.expect(P.rows == PP.rows and P.cols == PP.cols);

    for (0..A.rows) |i| for (0..A.cols) |j| {try std.testing.expect(@abs(A.at(i, j) - AA.at(i, j)) < 1e-14);};
    for (0..B.rows) |i| for (0..B.cols) |j| {try std.testing.expect(@abs(B.at(i, j) - BB.at(i, j)) < 1e-14);};
    for (0..C.rows) |i| for (0..C.cols) |j| {try std.testing.expect(@abs(C.at(i, j) - CC.at(i, j)) < 1e-14);};
    for (0..D.rows) |i| for (0..D.cols) |j| {try std.testing.expect(@abs(D.at(i, j) - DD.at(i, j)) < 1e-14);};
    for (0..E.rows) |i| for (0..E.cols) |j| {try std.testing.expect(@abs(E.at(i, j) - EE.at(i, j)) < 1e-14);};
    for (0..F.rows) |i| for (0..F.cols) |j| {try std.testing.expect(@abs(F.at(i, j) - FF.at(i, j)) < 1e-14);};
    for (0..G.rows) |i| for (0..G.cols) |j| {try std.testing.expect(@abs(G.at(i, j) - GG.at(i, j)) < 1e-14);};
    for (0..H.rows) |i| for (0..H.cols) |j| {try std.testing.expect(@abs(H.at(i, j) - HH.at(i, j)) < 1e-14);};
    for (0..I.rows) |i| for (0..I.cols) |j| {try std.testing.expect(@abs(I.at(i, j) - II.at(i, j)) < 1e-14);};
    for (0..J.rows) |i| for (0..J.cols) |j| {try std.testing.expect(@abs(J.at(i, j) - JJ.at(i, j)) < 1e-14);};
    for (0..K.rows) |i| for (0..K.cols) |j| {try std.testing.expect(@abs(K.at(i, j) - KK.at(i, j)) < 1e-14);};
    for (0..L.rows) |i| for (0..L.cols) |j| {try std.testing.expect(@abs(L.at(i, j) - LL.at(i, j)) < 1e-14);};
    for (0..M.rows) |i| for (0..M.cols) |j| {try std.testing.expect(@abs(M.at(i, j) - MM.at(i, j)) < 1e-14);};
    for (0..N.rows) |i| for (0..N.cols) |j| {try std.testing.expect(@abs(N.at(i, j) - NN.at(i, j)) < 1e-14);};
    for (0..O.rows) |i| for (0..O.cols) |j| {try std.testing.expect(@abs(O.at(i, j) - OO.at(i, j)) < 1e-14);};
    for (0..P.rows) |i| for (0..P.cols) |j| {try std.testing.expect(@abs(P.at(i, j) - PP.at(i, j)) < 1e-14);};
}

test "matrix_linspace" {
    const A = try matrix.Matrix(f64).init(1   , 1   , allocator); defer A.deinit();
    const B = try matrix.Matrix(f64).init(1   , 10  , allocator); defer B.deinit();
    const C = try matrix.Matrix(f64).init(1   , 100 , allocator); defer C.deinit();
    const D = try matrix.Matrix(f64).init(1   , 1000, allocator); defer D.deinit();
    const E = try matrix.Matrix(f64).init(10  , 1   , allocator); defer E.deinit();
    const F = try matrix.Matrix(f64).init(10  , 10  , allocator); defer F.deinit();
    const G = try matrix.Matrix(f64).init(10  , 100 , allocator); defer G.deinit();
    const H = try matrix.Matrix(f64).init(10  , 1000, allocator); defer H.deinit();
    const I = try matrix.Matrix(f64).init(100 , 1   , allocator); defer I.deinit();
    const J = try matrix.Matrix(f64).init(100 , 10  , allocator); defer J.deinit();
    const K = try matrix.Matrix(f64).init(100 , 100 , allocator); defer K.deinit();
    const L = try matrix.Matrix(f64).init(100 , 1000, allocator); defer L.deinit();
    const M = try matrix.Matrix(f64).init(1000, 1   , allocator); defer M.deinit();
    const N = try matrix.Matrix(f64).init(1000, 10  , allocator); defer N.deinit();
    const O = try matrix.Matrix(f64).init(1000, 100 , allocator); defer O.deinit();
    const P = try matrix.Matrix(f64).init(1000, 1000, allocator); defer P.deinit();

    A.linspace(0, helper.asfloat(f64, A.data.len) - 1); B.linspace(0, helper.asfloat(f64, B.data.len) - 1); C.linspace(0, helper.asfloat(f64, C.data.len) - 1); D.linspace(0, helper.asfloat(f64, D.data.len) - 1);
    E.linspace(0, helper.asfloat(f64, E.data.len) - 1); F.linspace(0, helper.asfloat(f64, F.data.len) - 1); G.linspace(0, helper.asfloat(f64, G.data.len) - 1); H.linspace(0, helper.asfloat(f64, H.data.len) - 1);
    I.linspace(0, helper.asfloat(f64, I.data.len) - 1); J.linspace(0, helper.asfloat(f64, J.data.len) - 1); K.linspace(0, helper.asfloat(f64, K.data.len) - 1); L.linspace(0, helper.asfloat(f64, L.data.len) - 1);
    M.linspace(0, helper.asfloat(f64, M.data.len) - 1); N.linspace(0, helper.asfloat(f64, N.data.len) - 1); O.linspace(0, helper.asfloat(f64, O.data.len) - 1); P.linspace(0, helper.asfloat(f64, P.data.len) - 1);

    for (0..A.rows) |i| for (0..A.cols) |j| {try std.testing.expect(A.at(i, j) == helper.asfloat(f64, i * A.cols + j));};
    for (0..B.rows) |i| for (0..B.cols) |j| {try std.testing.expect(B.at(i, j) == helper.asfloat(f64, i * B.cols + j));};
    for (0..C.rows) |i| for (0..C.cols) |j| {try std.testing.expect(C.at(i, j) == helper.asfloat(f64, i * C.cols + j));};
    for (0..D.rows) |i| for (0..D.cols) |j| {try std.testing.expect(D.at(i, j) == helper.asfloat(f64, i * D.cols + j));};
    for (0..E.rows) |i| for (0..E.cols) |j| {try std.testing.expect(E.at(i, j) == helper.asfloat(f64, i * E.cols + j));};
    for (0..F.rows) |i| for (0..F.cols) |j| {try std.testing.expect(F.at(i, j) == helper.asfloat(f64, i * F.cols + j));};
    for (0..G.rows) |i| for (0..G.cols) |j| {try std.testing.expect(G.at(i, j) == helper.asfloat(f64, i * G.cols + j));};
    for (0..H.rows) |i| for (0..H.cols) |j| {try std.testing.expect(H.at(i, j) == helper.asfloat(f64, i * H.cols + j));};
    for (0..I.rows) |i| for (0..I.cols) |j| {try std.testing.expect(I.at(i, j) == helper.asfloat(f64, i * I.cols + j));};
    for (0..J.rows) |i| for (0..J.cols) |j| {try std.testing.expect(J.at(i, j) == helper.asfloat(f64, i * J.cols + j));};
    for (0..K.rows) |i| for (0..K.cols) |j| {try std.testing.expect(K.at(i, j) == helper.asfloat(f64, i * K.cols + j));};
    for (0..L.rows) |i| for (0..L.cols) |j| {try std.testing.expect(L.at(i, j) == helper.asfloat(f64, i * L.cols + j));};
    for (0..M.rows) |i| for (0..M.cols) |j| {try std.testing.expect(M.at(i, j) == helper.asfloat(f64, i * M.cols + j));};
    for (0..N.rows) |i| for (0..N.cols) |j| {try std.testing.expect(N.at(i, j) == helper.asfloat(f64, i * N.cols + j));};
    for (0..O.rows) |i| for (0..O.cols) |j| {try std.testing.expect(O.at(i, j) == helper.asfloat(f64, i * O.cols + j));};
    for (0..P.rows) |i| for (0..P.cols) |j| {try std.testing.expect(P.at(i, j) == helper.asfloat(f64, i * P.cols + j));};
}

test "matrix_randn" {
    const A = try matrix.Matrix(f64).init(1   , 1   , allocator); defer A.deinit();
    const B = try matrix.Matrix(f64).init(1   , 10  , allocator); defer B.deinit();
    const C = try matrix.Matrix(f64).init(1   , 100 , allocator); defer C.deinit();
    const D = try matrix.Matrix(f64).init(1   , 1000, allocator); defer D.deinit();
    const E = try matrix.Matrix(f64).init(10  , 1   , allocator); defer E.deinit();
    const F = try matrix.Matrix(f64).init(10  , 10  , allocator); defer F.deinit();
    const G = try matrix.Matrix(f64).init(10  , 100 , allocator); defer G.deinit();
    const H = try matrix.Matrix(f64).init(10  , 1000, allocator); defer H.deinit();
    const I = try matrix.Matrix(f64).init(100 , 1   , allocator); defer I.deinit();
    const J = try matrix.Matrix(f64).init(100 , 10  , allocator); defer J.deinit();
    const K = try matrix.Matrix(f64).init(100 , 100 , allocator); defer K.deinit();
    const L = try matrix.Matrix(f64).init(100 , 1000, allocator); defer L.deinit();
    const M = try matrix.Matrix(f64).init(1000, 1   , allocator); defer M.deinit();
    const N = try matrix.Matrix(f64).init(1000, 10  , allocator); defer N.deinit();
    const O = try matrix.Matrix(f64).init(1000, 100 , allocator); defer O.deinit();
    const P = try matrix.Matrix(f64).init(1000, 1000, allocator); defer P.deinit();

    A.fill( 0      ); B.fill( 0      ); C.fill( 0       ); D.fill( 0       ); E.fill( 0       ); F.fill( 0       ); G.fill( 0       ); H.fill( 0       );
    I.fill( 0      ); J.fill( 0      ); K.fill( 0       ); L.fill( 0       ); M.fill( 0       ); N.fill( 0       ); O.fill( 0       ); P.fill( 0       );
    A.randn(0, 1, 0); B.randn(0, 1, 1); C.randn(0, 1, 2 ); D.randn(0, 1, 3 ); E.randn(0, 1, 4 ); F.randn(0, 1, 5 ); G.randn(0, 1, 6 ); H.randn(0, 1, 7 );
    I.randn(0, 1, 8); J.randn(0, 1, 9); K.randn(0, 1, 10); L.randn(0, 1, 11); M.randn(0, 1, 12); N.randn(0, 1, 13); O.randn(0, 1, 14); P.randn(0, 1, 15);

    for (0..A.rows) |i| for (0..A.cols) |j| {try std.testing.expect(A.at(i, j) != 0);};
    for (0..B.rows) |i| for (0..B.cols) |j| {try std.testing.expect(B.at(i, j) != 0);};
    for (0..C.rows) |i| for (0..C.cols) |j| {try std.testing.expect(C.at(i, j) != 0);};
    for (0..D.rows) |i| for (0..D.cols) |j| {try std.testing.expect(D.at(i, j) != 0);};
    for (0..E.rows) |i| for (0..E.cols) |j| {try std.testing.expect(E.at(i, j) != 0);};
    for (0..F.rows) |i| for (0..F.cols) |j| {try std.testing.expect(F.at(i, j) != 0);};
    for (0..G.rows) |i| for (0..G.cols) |j| {try std.testing.expect(G.at(i, j) != 0);};
    for (0..H.rows) |i| for (0..H.cols) |j| {try std.testing.expect(H.at(i, j) != 0);};
    for (0..I.rows) |i| for (0..I.cols) |j| {try std.testing.expect(I.at(i, j) != 0);};
    for (0..J.rows) |i| for (0..J.cols) |j| {try std.testing.expect(J.at(i, j) != 0);};
    for (0..K.rows) |i| for (0..K.cols) |j| {try std.testing.expect(K.at(i, j) != 0);};
    for (0..L.rows) |i| for (0..L.cols) |j| {try std.testing.expect(L.at(i, j) != 0);};
    for (0..M.rows) |i| for (0..M.cols) |j| {try std.testing.expect(M.at(i, j) != 0);};
    for (0..N.rows) |i| for (0..N.cols) |j| {try std.testing.expect(N.at(i, j) != 0);};
    for (0..O.rows) |i| for (0..O.cols) |j| {try std.testing.expect(O.at(i, j) != 0);};
    for (0..P.rows) |i| for (0..P.cols) |j| {try std.testing.expect(P.at(i, j) != 0);};
}

test "matrix_row" {
    const A = try matrix.Matrix(f64).init(1   , 1   , allocator); defer A.deinit();
    const B = try matrix.Matrix(f64).init(1   , 10  , allocator); defer B.deinit();
    const C = try matrix.Matrix(f64).init(1   , 100 , allocator); defer C.deinit();
    const D = try matrix.Matrix(f64).init(1   , 1000, allocator); defer D.deinit();
    const E = try matrix.Matrix(f64).init(10  , 1   , allocator); defer E.deinit();
    const F = try matrix.Matrix(f64).init(10  , 10  , allocator); defer F.deinit();
    const G = try matrix.Matrix(f64).init(10  , 100 , allocator); defer G.deinit();
    const H = try matrix.Matrix(f64).init(10  , 1000, allocator); defer H.deinit();
    const I = try matrix.Matrix(f64).init(100 , 1   , allocator); defer I.deinit();
    const J = try matrix.Matrix(f64).init(100 , 10  , allocator); defer J.deinit();
    const K = try matrix.Matrix(f64).init(100 , 100 , allocator); defer K.deinit();
    const L = try matrix.Matrix(f64).init(100 , 1000, allocator); defer L.deinit();
    const M = try matrix.Matrix(f64).init(1000, 1   , allocator); defer M.deinit();
    const N = try matrix.Matrix(f64).init(1000, 10  , allocator); defer N.deinit();
    const O = try matrix.Matrix(f64).init(1000, 100 , allocator); defer O.deinit();
    const P = try matrix.Matrix(f64).init(1000, 1000, allocator); defer P.deinit();

    A.randn(0, 1, 0); B.randn(0, 1, 1); C.randn(0, 1, 2 ); D.randn(0, 1, 3 ); E.randn(0, 1, 4 ); F.randn(0, 1, 5 ); G.randn(0, 1, 6 ); H.randn(0, 1, 7 );
    I.randn(0, 1, 8); J.randn(0, 1, 9); K.randn(0, 1, 10); L.randn(0, 1, 11); M.randn(0, 1, 12); N.randn(0, 1, 13); O.randn(0, 1, 14); P.randn(0, 1, 15);

    for (0..A.rows) |i| {const row = A.row(i); for (0..A.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == A.cols and row.at(0, j) == A.at(i, j));}
    for (0..B.rows) |i| {const row = B.row(i); for (0..B.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == B.cols and row.at(0, j) == B.at(i, j));}
    for (0..C.rows) |i| {const row = C.row(i); for (0..C.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == C.cols and row.at(0, j) == C.at(i, j));}
    for (0..D.rows) |i| {const row = D.row(i); for (0..D.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == D.cols and row.at(0, j) == D.at(i, j));}
    for (0..E.rows) |i| {const row = E.row(i); for (0..E.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == E.cols and row.at(0, j) == E.at(i, j));}
    for (0..F.rows) |i| {const row = F.row(i); for (0..F.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == F.cols and row.at(0, j) == F.at(i, j));}
    for (0..G.rows) |i| {const row = G.row(i); for (0..G.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == G.cols and row.at(0, j) == G.at(i, j));}
    for (0..H.rows) |i| {const row = H.row(i); for (0..H.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == H.cols and row.at(0, j) == H.at(i, j));}
    for (0..I.rows) |i| {const row = I.row(i); for (0..I.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == I.cols and row.at(0, j) == I.at(i, j));}
    for (0..J.rows) |i| {const row = J.row(i); for (0..J.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == J.cols and row.at(0, j) == J.at(i, j));}
    for (0..K.rows) |i| {const row = K.row(i); for (0..K.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == K.cols and row.at(0, j) == K.at(i, j));}
    for (0..L.rows) |i| {const row = L.row(i); for (0..L.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == L.cols and row.at(0, j) == L.at(i, j));}
    for (0..M.rows) |i| {const row = M.row(i); for (0..M.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == M.cols and row.at(0, j) == M.at(i, j));}
    for (0..N.rows) |i| {const row = N.row(i); for (0..N.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == N.cols and row.at(0, j) == N.at(i, j));}
    for (0..O.rows) |i| {const row = O.row(i); for (0..O.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == O.cols and row.at(0, j) == O.at(i, j));}
    for (0..P.rows) |i| {const row = P.row(i); for (0..P.cols) |j| try std.testing.expect(row.rows == 1 and row.cols == P.cols and row.at(0, j) == P.at(i, j));}
}

test "matrix_sa" {
    const A = try matrix.Matrix(f64).init(1   , 1   , allocator); defer A.deinit();
    const B = try matrix.Matrix(f64).init(1   , 10  , allocator); defer B.deinit();
    const C = try matrix.Matrix(f64).init(1   , 100 , allocator); defer C.deinit();
    const D = try matrix.Matrix(f64).init(1   , 1000, allocator); defer D.deinit();
    const E = try matrix.Matrix(f64).init(10  , 1   , allocator); defer E.deinit();
    const F = try matrix.Matrix(f64).init(10  , 10  , allocator); defer F.deinit();
    const G = try matrix.Matrix(f64).init(10  , 100 , allocator); defer G.deinit();
    const H = try matrix.Matrix(f64).init(10  , 1000, allocator); defer H.deinit();
    const I = try matrix.Matrix(f64).init(100 , 1   , allocator); defer I.deinit();
    const J = try matrix.Matrix(f64).init(100 , 10  , allocator); defer J.deinit();
    const K = try matrix.Matrix(f64).init(100 , 100 , allocator); defer K.deinit();
    const L = try matrix.Matrix(f64).init(100 , 1000, allocator); defer L.deinit();
    const M = try matrix.Matrix(f64).init(1000, 1   , allocator); defer M.deinit();
    const N = try matrix.Matrix(f64).init(1000, 10  , allocator); defer N.deinit();
    const O = try matrix.Matrix(f64).init(1000, 100 , allocator); defer O.deinit();
    const P = try matrix.Matrix(f64).init(1000, 1000, allocator); defer P.deinit();

    A.randn(0, 1, 0); B.randn(0, 1, 1); C.randn(0, 1, 2 ); D.randn(0, 1, 3 ); E.randn(0, 1, 4 ); F.randn(0, 1, 5 ); G.randn(0, 1, 6 ); H.randn(0, 1, 7 );
    I.randn(0, 1, 8); J.randn(0, 1, 9); K.randn(0, 1, 10); L.randn(0, 1, 11); M.randn(0, 1, 12); N.randn(0, 1, 13); O.randn(0, 1, 14); P.randn(0, 1, 15);

    const AA = A.sa();
    const BB = B.sa();
    const CC = C.sa();
    const DD = D.sa();
    const EE = E.sa();
    const FF = F.sa();
    const GG = G.sa();
    const HH = H.sa();
    const II = I.sa();
    const JJ = J.sa();
    const KK = K.sa();
    const LL = L.sa();
    const MM = M.sa();
    const NN = N.sa();
    const OO = O.sa();
    const PP = P.sa();

    for (0..A.rows) |i| for (0..A.cols) |j| {try std.testing.expect(A.at(i, j) == AA.at(i * A.cols + j));};
    for (0..B.rows) |i| for (0..B.cols) |j| {try std.testing.expect(B.at(i, j) == BB.at(i * B.cols + j));};
    for (0..C.rows) |i| for (0..C.cols) |j| {try std.testing.expect(C.at(i, j) == CC.at(i * C.cols + j));};
    for (0..D.rows) |i| for (0..D.cols) |j| {try std.testing.expect(D.at(i, j) == DD.at(i * D.cols + j));};
    for (0..E.rows) |i| for (0..E.cols) |j| {try std.testing.expect(E.at(i, j) == EE.at(i * E.cols + j));};
    for (0..F.rows) |i| for (0..F.cols) |j| {try std.testing.expect(F.at(i, j) == FF.at(i * F.cols + j));};
    for (0..G.rows) |i| for (0..G.cols) |j| {try std.testing.expect(G.at(i, j) == GG.at(i * G.cols + j));};
    for (0..H.rows) |i| for (0..H.cols) |j| {try std.testing.expect(H.at(i, j) == HH.at(i * H.cols + j));};
    for (0..I.rows) |i| for (0..I.cols) |j| {try std.testing.expect(I.at(i, j) == II.at(i * I.cols + j));};
    for (0..J.rows) |i| for (0..J.cols) |j| {try std.testing.expect(J.at(i, j) == JJ.at(i * J.cols + j));};
    for (0..K.rows) |i| for (0..K.cols) |j| {try std.testing.expect(K.at(i, j) == KK.at(i * K.cols + j));};
    for (0..L.rows) |i| for (0..L.cols) |j| {try std.testing.expect(L.at(i, j) == LL.at(i * L.cols + j));};
    for (0..M.rows) |i| for (0..M.cols) |j| {try std.testing.expect(M.at(i, j) == MM.at(i * M.cols + j));};
    for (0..N.rows) |i| for (0..N.cols) |j| {try std.testing.expect(N.at(i, j) == NN.at(i * N.cols + j));};
    for (0..O.rows) |i| for (0..O.cols) |j| {try std.testing.expect(O.at(i, j) == OO.at(i * O.cols + j));};
    for (0..P.rows) |i| for (0..P.cols) |j| {try std.testing.expect(P.at(i, j) == PP.at(i * P.cols + j));};
}

test "matrix_vector" {
    const A = try matrix.Matrix(f64).init(1   , 1   , allocator); defer A.deinit();
    const B = try matrix.Matrix(f64).init(1   , 10  , allocator); defer B.deinit();
    const C = try matrix.Matrix(f64).init(1   , 100 , allocator); defer C.deinit();
    const D = try matrix.Matrix(f64).init(1   , 1000, allocator); defer D.deinit();
    const E = try matrix.Matrix(f64).init(10  , 1   , allocator); defer E.deinit();
    const F = try matrix.Matrix(f64).init(10  , 10  , allocator); defer F.deinit();
    const G = try matrix.Matrix(f64).init(10  , 100 , allocator); defer G.deinit();
    const H = try matrix.Matrix(f64).init(10  , 1000, allocator); defer H.deinit();
    const I = try matrix.Matrix(f64).init(100 , 1   , allocator); defer I.deinit();
    const J = try matrix.Matrix(f64).init(100 , 10  , allocator); defer J.deinit();
    const K = try matrix.Matrix(f64).init(100 , 100 , allocator); defer K.deinit();
    const L = try matrix.Matrix(f64).init(100 , 1000, allocator); defer L.deinit();
    const M = try matrix.Matrix(f64).init(1000, 1   , allocator); defer M.deinit();
    const N = try matrix.Matrix(f64).init(1000, 10  , allocator); defer N.deinit();
    const O = try matrix.Matrix(f64).init(1000, 100 , allocator); defer O.deinit();
    const P = try matrix.Matrix(f64).init(1000, 1000, allocator); defer P.deinit();

    A.randn(0, 1, 0); B.randn(0, 1, 1); C.randn(0, 1, 2 ); D.randn(0, 1, 3 ); E.randn(0, 1, 4 ); F.randn(0, 1, 5 ); G.randn(0, 1, 6 ); H.randn(0, 1, 7 );
    I.randn(0, 1, 8); J.randn(0, 1, 9); K.randn(0, 1, 10); L.randn(0, 1, 11); M.randn(0, 1, 12); N.randn(0, 1, 13); O.randn(0, 1, 14); P.randn(0, 1, 15);

    const AA = A.vector();
    const BB = B.vector();
    const CC = C.vector();
    const DD = D.vector();
    const EE = E.vector();
    const FF = F.vector();
    const GG = G.vector();
    const HH = H.vector();
    const II = I.vector();
    const JJ = J.vector();
    const KK = K.vector();
    const LL = L.vector();
    const MM = M.vector();
    const NN = N.vector();
    const OO = O.vector();
    const PP = P.vector();

    for (0..A.rows) |i| for (0..A.cols) |j| {try std.testing.expect(A.at(i, j) == AA.at(i * A.cols + j));};
    for (0..B.rows) |i| for (0..B.cols) |j| {try std.testing.expect(B.at(i, j) == BB.at(i * B.cols + j));};
    for (0..C.rows) |i| for (0..C.cols) |j| {try std.testing.expect(C.at(i, j) == CC.at(i * C.cols + j));};
    for (0..D.rows) |i| for (0..D.cols) |j| {try std.testing.expect(D.at(i, j) == DD.at(i * D.cols + j));};
    for (0..E.rows) |i| for (0..E.cols) |j| {try std.testing.expect(E.at(i, j) == EE.at(i * E.cols + j));};
    for (0..F.rows) |i| for (0..F.cols) |j| {try std.testing.expect(F.at(i, j) == FF.at(i * F.cols + j));};
    for (0..G.rows) |i| for (0..G.cols) |j| {try std.testing.expect(G.at(i, j) == GG.at(i * G.cols + j));};
    for (0..H.rows) |i| for (0..H.cols) |j| {try std.testing.expect(H.at(i, j) == HH.at(i * H.cols + j));};
    for (0..I.rows) |i| for (0..I.cols) |j| {try std.testing.expect(I.at(i, j) == II.at(i * I.cols + j));};
    for (0..J.rows) |i| for (0..J.cols) |j| {try std.testing.expect(J.at(i, j) == JJ.at(i * J.cols + j));};
    for (0..K.rows) |i| for (0..K.cols) |j| {try std.testing.expect(K.at(i, j) == KK.at(i * K.cols + j));};
    for (0..L.rows) |i| for (0..L.cols) |j| {try std.testing.expect(L.at(i, j) == LL.at(i * L.cols + j));};
    for (0..M.rows) |i| for (0..M.cols) |j| {try std.testing.expect(M.at(i, j) == MM.at(i * M.cols + j));};
    for (0..N.rows) |i| for (0..N.cols) |j| {try std.testing.expect(N.at(i, j) == NN.at(i * N.cols + j));};
    for (0..O.rows) |i| for (0..O.cols) |j| {try std.testing.expect(O.at(i, j) == OO.at(i * O.cols + j));};
    for (0..P.rows) |i| for (0..P.cols) |j| {try std.testing.expect(P.at(i, j) == PP.at(i * P.cols + j));};
}
