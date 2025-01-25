const std = @import("std");

const helper = @import("acorn").helper;
const vector = @import("acorn").vector;

const Vector = @import("acorn").Vector;

const allocator = std.testing.allocator;

test "vector_at" {
    const u = try vector.Vector(f64).init(1      , allocator); defer u.deinit();
    const v = try vector.Vector(f64).init(10     , allocator); defer v.deinit();
    const w = try vector.Vector(f64).init(100    , allocator); defer w.deinit();
    const x = try vector.Vector(f64).init(1000   , allocator); defer x.deinit();
    const y = try vector.Vector(f64).init(10000  , allocator); defer y.deinit();
    const z = try vector.Vector(f64).init(100000 , allocator); defer z.deinit();
    const a = try vector.Vector(f64).init(1000000, allocator); defer a.deinit();

    for (0..u.rows) |i| u.ptr(i).* = helper.asfloat(f64, i);
    for (0..v.rows) |i| v.ptr(i).* = helper.asfloat(f64, i);
    for (0..w.rows) |i| w.ptr(i).* = helper.asfloat(f64, i);
    for (0..x.rows) |i| x.ptr(i).* = helper.asfloat(f64, i);
    for (0..y.rows) |i| y.ptr(i).* = helper.asfloat(f64, i);
    for (0..z.rows) |i| z.ptr(i).* = helper.asfloat(f64, i);
    for (0..a.rows) |i| a.ptr(i).* = helper.asfloat(f64, i);

    for (0..u.rows) |i| try std.testing.expect(u.at(i) == helper.asfloat(f64, i));
    for (0..v.rows) |i| try std.testing.expect(v.at(i) == helper.asfloat(f64, i));
    for (0..w.rows) |i| try std.testing.expect(w.at(i) == helper.asfloat(f64, i));
    for (0..x.rows) |i| try std.testing.expect(x.at(i) == helper.asfloat(f64, i));
    for (0..y.rows) |i| try std.testing.expect(y.at(i) == helper.asfloat(f64, i));
    for (0..z.rows) |i| try std.testing.expect(z.at(i) == helper.asfloat(f64, i));
    for (0..a.rows) |i| try std.testing.expect(a.at(i) == helper.asfloat(f64, i));
}

test "vector_clone" {
    const u = try vector.Vector(f64).init(1      , allocator); defer u.deinit();
    const v = try vector.Vector(f64).init(10     , allocator); defer v.deinit();
    const w = try vector.Vector(f64).init(100    , allocator); defer w.deinit();
    const x = try vector.Vector(f64).init(1000   , allocator); defer x.deinit();
    const y = try vector.Vector(f64).init(10000  , allocator); defer y.deinit();
    const z = try vector.Vector(f64).init(100000 , allocator); defer z.deinit();
    const a = try vector.Vector(f64).init(1000000, allocator); defer a.deinit();

    u.randn(0, 1, 0); v.randn(0, 1, 1); w.randn(0, 1, 2); x.randn(0, 1, 3); y.randn(0, 1, 4); z.randn(0, 1, 5); a.randn(0, 1, 6);

    const uu = try u.clone(); defer uu.deinit();
    const vv = try v.clone(); defer vv.deinit();
    const ww = try w.clone(); defer ww.deinit();
    const xx = try x.clone(); defer xx.deinit();
    const yy = try y.clone(); defer yy.deinit();
    const zz = try z.clone(); defer zz.deinit();
    const aa = try a.clone(); defer aa.deinit();

    for (0..u.rows) |i| try std.testing.expect(u.at(i) == uu.at(i));
    for (0..v.rows) |i| try std.testing.expect(v.at(i) == vv.at(i));
    for (0..w.rows) |i| try std.testing.expect(w.at(i) == ww.at(i));
    for (0..x.rows) |i| try std.testing.expect(x.at(i) == xx.at(i));
    for (0..y.rows) |i| try std.testing.expect(y.at(i) == yy.at(i));
    for (0..z.rows) |i| try std.testing.expect(z.at(i) == zz.at(i));
    for (0..a.rows) |i| try std.testing.expect(a.at(i) == aa.at(i));
}

test "vector_complex" {
    const u = try vector.Vector(f64).init(1      , allocator); defer u.deinit();
    const v = try vector.Vector(f64).init(10     , allocator); defer v.deinit();
    const w = try vector.Vector(f64).init(100    , allocator); defer w.deinit();
    const x = try vector.Vector(f64).init(1000   , allocator); defer x.deinit();
    const y = try vector.Vector(f64).init(10000  , allocator); defer y.deinit();
    const z = try vector.Vector(f64).init(100000 , allocator); defer z.deinit();
    const a = try vector.Vector(f64).init(1000000, allocator); defer a.deinit();

    u.randn(0, 1, 0); v.randn(0, 1, 1); w.randn(0, 1, 2); x.randn(0, 1, 3); y.randn(0, 1, 4); z.randn(0, 1, 5); a.randn(0, 1, 6);

    const uu = try u.complex(); defer uu.deinit();
    const vv = try v.complex(); defer vv.deinit();
    const ww = try w.complex(); defer ww.deinit();
    const xx = try x.complex(); defer xx.deinit();
    const yy = try y.complex(); defer yy.deinit();
    const zz = try z.complex(); defer zz.deinit();
    const aa = try a.complex(); defer aa.deinit();

    for (0..u.rows) |i| try std.testing.expect(u.at(i) == uu.at(i).re and uu.at(i).im == 0);
    for (0..v.rows) |i| try std.testing.expect(v.at(i) == vv.at(i).re and vv.at(i).im == 0);
    for (0..w.rows) |i| try std.testing.expect(w.at(i) == ww.at(i).re and ww.at(i).im == 0);
    for (0..x.rows) |i| try std.testing.expect(x.at(i) == xx.at(i).re and xx.at(i).im == 0);
    for (0..y.rows) |i| try std.testing.expect(y.at(i) == yy.at(i).re and yy.at(i).im == 0);
    for (0..z.rows) |i| try std.testing.expect(z.at(i) == zz.at(i).re and zz.at(i).im == 0);
    for (0..a.rows) |i| try std.testing.expect(a.at(i) == aa.at(i).re and aa.at(i).im == 0);
}

test "vector_fill" {
    const u = try vector.Vector(f64).init(1      , allocator); defer u.deinit();
    const v = try vector.Vector(f64).init(10     , allocator); defer v.deinit();
    const w = try vector.Vector(f64).init(100    , allocator); defer w.deinit();
    const x = try vector.Vector(f64).init(1000   , allocator); defer x.deinit();
    const y = try vector.Vector(f64).init(10000  , allocator); defer y.deinit();
    const z = try vector.Vector(f64).init(100000 , allocator); defer z.deinit();
    const a = try vector.Vector(f64).init(1000000, allocator); defer a.deinit();

    u.fill(0); v.fill(1); w.fill(2); x.fill(3); y.fill(4); z.fill(5); a.fill(6);

    for (0..u.rows) |i| try std.testing.expect(u.at(i) == 0);
    for (0..v.rows) |i| try std.testing.expect(v.at(i) == 1);
    for (0..w.rows) |i| try std.testing.expect(w.at(i) == 2);
    for (0..x.rows) |i| try std.testing.expect(x.at(i) == 3);
    for (0..y.rows) |i| try std.testing.expect(y.at(i) == 4);
    for (0..z.rows) |i| try std.testing.expect(z.at(i) == 5);
    for (0..a.rows) |i| try std.testing.expect(a.at(i) == 6);
}

test "vector_init" {
    const u = try vector.Vector(f64).init(1      , allocator); defer u.deinit();
    const v = try vector.Vector(f64).init(10     , allocator); defer v.deinit();
    const w = try vector.Vector(f64).init(100    , allocator); defer w.deinit();
    const x = try vector.Vector(f64).init(1000   , allocator); defer x.deinit();
    const y = try vector.Vector(f64).init(10000  , allocator); defer y.deinit();
    const z = try vector.Vector(f64).init(100000 , allocator); defer z.deinit();
    const a = try vector.Vector(f64).init(1000000, allocator); defer a.deinit();

    try std.testing.expect(u.rows == 1      );
    try std.testing.expect(v.rows == 10     );
    try std.testing.expect(w.rows == 100    );
    try std.testing.expect(x.rows == 1000   );
    try std.testing.expect(y.rows == 10000  );
    try std.testing.expect(z.rows == 100000 );
    try std.testing.expect(a.rows == 1000000);
}

test "vector_matrix" {
    const u = try vector.Vector(f64).init(1      , allocator); defer u.deinit();
    const v = try vector.Vector(f64).init(10     , allocator); defer v.deinit();
    const w = try vector.Vector(f64).init(100    , allocator); defer w.deinit();
    const x = try vector.Vector(f64).init(1000   , allocator); defer x.deinit();
    const y = try vector.Vector(f64).init(10000  , allocator); defer y.deinit();
    const z = try vector.Vector(f64).init(100000 , allocator); defer z.deinit();
    const a = try vector.Vector(f64).init(1000000, allocator); defer a.deinit();

    u.randn(0, 1, 0); v.randn(0, 1, 1); w.randn(0, 1, 2); x.randn(0, 1, 3); y.randn(0, 1, 4); z.randn(0, 1, 5); a.randn(0, 1, 6);

    const uu = u.matrix();
    const vv = v.matrix();
    const ww = w.matrix();
    const xx = x.matrix();
    const yy = y.matrix();
    const zz = z.matrix();
    const aa = a.matrix();

    try std.testing.expect(uu.rows == u.rows and uu.cols == 1);
    try std.testing.expect(vv.rows == v.rows and vv.cols == 1);
    try std.testing.expect(ww.rows == w.rows and ww.cols == 1);
    try std.testing.expect(xx.rows == x.rows and xx.cols == 1);
    try std.testing.expect(yy.rows == y.rows and yy.cols == 1);
    try std.testing.expect(zz.rows == z.rows and zz.cols == 1);
    try std.testing.expect(aa.rows == a.rows and aa.cols == 1);

    for (0..u.rows) |i| try std.testing.expect(u.at(i) == uu.at(i, 0));
    for (0..v.rows) |i| try std.testing.expect(v.at(i) == vv.at(i, 0));
    for (0..w.rows) |i| try std.testing.expect(w.at(i) == ww.at(i, 0));
    for (0..x.rows) |i| try std.testing.expect(x.at(i) == xx.at(i, 0));
    for (0..y.rows) |i| try std.testing.expect(y.at(i) == yy.at(i, 0));
    for (0..z.rows) |i| try std.testing.expect(z.at(i) == zz.at(i, 0));
    for (0..a.rows) |i| try std.testing.expect(a.at(i) == aa.at(i, 0));
}

test "vector_randn" {
    const u = try vector.Vector(f64).init(1      , allocator); defer u.deinit();
    const v = try vector.Vector(f64).init(10     , allocator); defer v.deinit();
    const w = try vector.Vector(f64).init(100    , allocator); defer w.deinit();
    const x = try vector.Vector(f64).init(1000   , allocator); defer x.deinit();
    const y = try vector.Vector(f64).init(10000  , allocator); defer y.deinit();
    const z = try vector.Vector(f64).init(100000 , allocator); defer z.deinit();
    const a = try vector.Vector(f64).init(1000000, allocator); defer a.deinit();

    u.fill( 0      ); v.fill( 0      ); w.fill( 0      ); x.fill( 0      ); y.fill( 0      ); z.fill( 0      ); a.fill( 0      );
    u.randn(0, 1, 0); v.randn(0, 1, 1); w.randn(0, 1, 2); x.randn(0, 1, 3); y.randn(0, 1, 4); z.randn(0, 1, 5); a.randn(0, 1, 6);

    for (0..u.rows) |i| try std.testing.expect(u.at(i) != 0);
    for (0..v.rows) |i| try std.testing.expect(v.at(i) != 0);
    for (0..w.rows) |i| try std.testing.expect(w.at(i) != 0);
    for (0..x.rows) |i| try std.testing.expect(x.at(i) != 0);
    for (0..y.rows) |i| try std.testing.expect(y.at(i) != 0);
    for (0..z.rows) |i| try std.testing.expect(z.at(i) != 0);
    for (0..a.rows) |i| try std.testing.expect(a.at(i) != 0);
}
