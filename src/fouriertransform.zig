const std = @import("std"); const Complex = std.math.Complex;

const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;
const bitrev  = @import("helper.zig").bitrev ;

pub fn fft(comptime T: type, out: []Complex(T), in: []Complex(T), factor: i32) void {
    const n = in.len; const logn: u6 = @intCast(std.math.log2(n)); for (0..n) |i| out[bitrev(i, logn)] = in[i];

    for (0..logn) |i| {

        const m = std.math.pow(usize, 2, i + 1); const mix = std.math.complex.exp(Complex(T).init(0, asfloat(T, factor) * 2 * std.math.pi / asfloat(T, m)));

        for (0..n / m) |j| {

            var omega = Complex(T).init(1, 0);

            for (0..m / 2) |k| {

                const t = omega.mul(out[j * m + k + m / 2]); const u = out[j * m + k];

                out[j * m + k        ] = u.add(t);
                out[j * m + k + m / 2] = u.sub(t);

                omega = omega.mul(mix);
            }
        }
    }

    if (factor > 0) for (0..out.len) |i| {out[i] = out[i].div(Complex(T).init(asfloat(T, out.len), 0));};
}

test "fft" {
    var A = try Vector(Complex(f64)).init(4, std.testing.allocator); defer A.deinit();
    var B = try Vector(Complex(f64)).init(4, std.testing.allocator); defer B.deinit();

    A.ptr(0).* = Complex(f64).init(3, 0);
    A.ptr(1).* = Complex(f64).init(6, 0);
    A.ptr(2).* = Complex(f64).init(1, 0);
    A.ptr(3).* = Complex(f64).init(8, 0);

    fft(f64, B.data, A.data, -1);

    for (0..A.data.len) |i| {
        std.debug.print("{d}: {d:12.6} {d:12.6}\n", .{i + 1, A.data[i].re, A.data[i].im});
    }

    std.debug.print("\n", .{});

    for (0..B.data.len) |i| {
        std.debug.print("{d}: {d:12.6} {d:12.6}\n", .{i + 1, B.data[i].re, B.data[i].im});
    }

    try std.testing.expect(true);
}
