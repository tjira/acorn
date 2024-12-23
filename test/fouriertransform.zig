const std = @import("std"); const Complex = std.math.complex.Complex;

const ftr = @import("acorn").ftr;
const vec = @import("acorn").vec;

const Vector = @import("acorn").Vector;

test "ftr_fft" {
    var u: Vector(Complex(f64)) = undefined;
    var v: Vector(Complex(f64)) = undefined;
    var w: Vector(Complex(f64)) = undefined;

    u = try @TypeOf(u).init(1, std.testing.allocator);
    v = try @TypeOf(v).init(1, std.testing.allocator);
    w = try @TypeOf(w).init(1, std.testing.allocator);

    u.ptr(0).* = @TypeOf(u.at(0)).init(1, 0);
    w.ptr(0).* = @TypeOf(w.at(0)).init(1, 0);

    ftr.fft(@TypeOf(u.at(0).re), v.sa(), u.sa(), -1); try std.testing.expect(vec.eq(@TypeOf(u.at(0)), v, w, 1e-12));
    ftr.fft(@TypeOf(u.at(0).re), w.sa(), v.sa(),  1); try std.testing.expect(vec.eq(@TypeOf(u.at(0)), w, u, 1e-12));

    u.deinit(); v.deinit(); w.deinit();

    u = try @TypeOf(u).init(2, std.testing.allocator);
    v = try @TypeOf(v).init(2, std.testing.allocator);
    w = try @TypeOf(w).init(2, std.testing.allocator);

    u.ptr(0).* = @TypeOf(u.at(0)).init( 1, 0);
    u.ptr(1).* = @TypeOf(u.at(0)).init( 2, 0);
    w.ptr(0).* = @TypeOf(w.at(0)).init( 3, 0);
    w.ptr(1).* = @TypeOf(w.at(0)).init(-1, 0);

    ftr.fft(@TypeOf(u.at(0).re), v.sa(), u.sa(), -1); try std.testing.expect(vec.eq(@TypeOf(u.at(0)), v, w, 1e-12));
    ftr.fft(@TypeOf(u.at(0).re), w.sa(), v.sa(),  1); try std.testing.expect(vec.eq(@TypeOf(u.at(0)), w, u, 1e-12));

    u.deinit(); v.deinit(); w.deinit();

    u = try @TypeOf(u).init(4, std.testing.allocator);
    v = try @TypeOf(v).init(4, std.testing.allocator);
    w = try @TypeOf(w).init(4, std.testing.allocator);

    u.ptr(0).* = @TypeOf(u.at(0)).init( 1,  0);
    u.ptr(1).* = @TypeOf(u.at(0)).init( 2,  0);
    u.ptr(2).* = @TypeOf(u.at(0)).init( 3,  0);
    u.ptr(3).* = @TypeOf(u.at(0)).init( 4,  0);
    w.ptr(0).* = @TypeOf(w.at(0)).init(10,  0);
    w.ptr(1).* = @TypeOf(w.at(0)).init(-2,  2);
    w.ptr(2).* = @TypeOf(w.at(0)).init(-2,  0);
    w.ptr(3).* = @TypeOf(w.at(0)).init(-2, -2);

    ftr.fft(@TypeOf(u.at(0).re), v.sa(), u.sa(), -1); try std.testing.expect(vec.eq(@TypeOf(u.at(0)), v, w, 1e-12));
    ftr.fft(@TypeOf(u.at(0).re), w.sa(), v.sa(),  1); try std.testing.expect(vec.eq(@TypeOf(u.at(0)), w, u, 1e-12));

    u.deinit(); v.deinit(); w.deinit();
}
