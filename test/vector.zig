const std = @import("std");

const vec = @import("acorn").vec;

const Vector = @import("acorn").Vector;

const allocator = std.testing.allocator;

test "vec_init" {
    const u = try vec.Vector(f64).init(1      , allocator); defer u.deinit();
    const v = try vec.Vector(f64).init(10     , allocator); defer v.deinit();
    const w = try vec.Vector(f64).init(100    , allocator); defer w.deinit();
    const x = try vec.Vector(f64).init(1000   , allocator); defer x.deinit();
    const y = try vec.Vector(f64).init(10000  , allocator); defer y.deinit();
    const z = try vec.Vector(f64).init(100000 , allocator); defer z.deinit();
    const a = try vec.Vector(f64).init(1000000, allocator); defer a.deinit();

    try std.testing.expect(u.rows == 1      );
    try std.testing.expect(v.rows == 10     );
    try std.testing.expect(w.rows == 100    );
    try std.testing.expect(x.rows == 1000   );
    try std.testing.expect(y.rows == 10000  );
    try std.testing.expect(z.rows == 100000 );
    try std.testing.expect(a.rows == 1000000);
}
