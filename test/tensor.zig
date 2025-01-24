const std = @import("std");

const tensor = @import("acorn").tensor;

const Tensor = @import("acorn").Tensor;

const allocator = std.testing.allocator;

test "ten_init" {
    const A = try tensor.Tensor(f64).init(&[_]usize{1         }, allocator); defer A.deinit();
    const B = try tensor.Tensor(f64).init(&[_]usize{10        }, allocator); defer B.deinit();
    const C = try tensor.Tensor(f64).init(&[_]usize{100       }, allocator); defer C.deinit();
    const D = try tensor.Tensor(f64).init(&[_]usize{1000      }, allocator); defer D.deinit();
    const E = try tensor.Tensor(f64).init(&[_]usize{10000     }, allocator); defer E.deinit();
    const F = try tensor.Tensor(f64).init(&[_]usize{100000    }, allocator); defer F.deinit();
    const G = try tensor.Tensor(f64).init(&[_]usize{1000000   }, allocator); defer G.deinit();
    const H = try tensor.Tensor(f64).init(&[_]usize{1, 1      }, allocator); defer H.deinit();
    const I = try tensor.Tensor(f64).init(&[_]usize{1, 10     }, allocator); defer I.deinit();
    const J = try tensor.Tensor(f64).init(&[_]usize{1, 100    }, allocator); defer J.deinit();
    const K = try tensor.Tensor(f64).init(&[_]usize{1, 1000   }, allocator); defer K.deinit();
    const L = try tensor.Tensor(f64).init(&[_]usize{10, 1     }, allocator); defer L.deinit();
    const M = try tensor.Tensor(f64).init(&[_]usize{10, 10    }, allocator); defer M.deinit();
    const N = try tensor.Tensor(f64).init(&[_]usize{10, 100   }, allocator); defer N.deinit();
    const O = try tensor.Tensor(f64).init(&[_]usize{10, 1000  }, allocator); defer O.deinit();
    const P = try tensor.Tensor(f64).init(&[_]usize{100, 1    }, allocator); defer P.deinit();
    const Q = try tensor.Tensor(f64).init(&[_]usize{100, 10   }, allocator); defer Q.deinit();
    const R = try tensor.Tensor(f64).init(&[_]usize{100, 100  }, allocator); defer R.deinit();
    const S = try tensor.Tensor(f64).init(&[_]usize{100, 1000 }, allocator); defer S.deinit();
    const T = try tensor.Tensor(f64).init(&[_]usize{1000, 1   }, allocator); defer T.deinit();
    const U = try tensor.Tensor(f64).init(&[_]usize{1000, 10  }, allocator); defer U.deinit();
    const V = try tensor.Tensor(f64).init(&[_]usize{1000, 100 }, allocator); defer V.deinit();
    const W = try tensor.Tensor(f64).init(&[_]usize{1000, 1000}, allocator); defer W.deinit();

    try std.testing.expect(A.shape[0] == 1                             );
    try std.testing.expect(B.shape[0] == 10                            );
    try std.testing.expect(C.shape[0] == 100                           );
    try std.testing.expect(D.shape[0] == 1000                          );
    try std.testing.expect(E.shape[0] == 10000                         );
    try std.testing.expect(F.shape[0] == 100000                        );
    try std.testing.expect(G.shape[0] == 1000000                       );
    try std.testing.expect(H.shape[0] == 1       and H.shape[1] == 1   );
    try std.testing.expect(I.shape[0] == 1       and I.shape[1] == 10  );
    try std.testing.expect(J.shape[0] == 1       and J.shape[1] == 100 );
    try std.testing.expect(K.shape[0] == 1       and K.shape[1] == 1000);
    try std.testing.expect(L.shape[0] == 10      and L.shape[1] == 1   );
    try std.testing.expect(M.shape[0] == 10      and M.shape[1] == 10  );
    try std.testing.expect(N.shape[0] == 10      and N.shape[1] == 100 );
    try std.testing.expect(O.shape[0] == 10      and O.shape[1] == 1000);
    try std.testing.expect(P.shape[0] == 100     and P.shape[1] == 1   );
    try std.testing.expect(Q.shape[0] == 100     and Q.shape[1] == 10  );
    try std.testing.expect(R.shape[0] == 100     and R.shape[1] == 100 );
    try std.testing.expect(S.shape[0] == 100     and S.shape[1] == 1000);
    try std.testing.expect(T.shape[0] == 1000    and T.shape[1] == 1   );
    try std.testing.expect(U.shape[0] == 1000    and U.shape[1] == 10  );
    try std.testing.expect(V.shape[0] == 1000    and V.shape[1] == 100 );
    try std.testing.expect(W.shape[0] == 1000    and W.shape[1] == 1000);
}
