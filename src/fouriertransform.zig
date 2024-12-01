const std = @import("std"); const Complex = std.math.Complex;

const vec = @import("vector.zig");

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
