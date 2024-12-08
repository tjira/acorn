const std = @import("std"); const Complex = std.math.complex.Complex;

const vec = @import("acorn").vec;

const Vector = @import("acorn").Vector;

test "vec_init" {
    var u: Vector(f64) = undefined; var v: Vector(Complex(f64)) = undefined;

    u = try @TypeOf(u).init(1, std.testing.allocator); try std.testing.expect(u.rows == 1 and u.data.len == 1); u.deinit();
    u = try @TypeOf(u).init(2, std.testing.allocator); try std.testing.expect(u.rows == 2 and u.data.len == 2); u.deinit();
    u = try @TypeOf(u).init(3, std.testing.allocator); try std.testing.expect(u.rows == 3 and u.data.len == 3); u.deinit();
    u = try @TypeOf(u).init(4, std.testing.allocator); try std.testing.expect(u.rows == 4 and u.data.len == 4); u.deinit();

    v = try @TypeOf(v).init(1, std.testing.allocator); try std.testing.expect(v.rows == 1 and v.data.len == 1); v.deinit();
    v = try @TypeOf(v).init(2, std.testing.allocator); try std.testing.expect(v.rows == 2 and v.data.len == 2); v.deinit();
    v = try @TypeOf(v).init(3, std.testing.allocator); try std.testing.expect(v.rows == 3 and v.data.len == 3); v.deinit();
    v = try @TypeOf(v).init(4, std.testing.allocator); try std.testing.expect(v.rows == 4 and v.data.len == 4); v.deinit();
}
