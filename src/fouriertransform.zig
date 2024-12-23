const std = @import("std"); const Complex = std.math.Complex;

const vec = @import("vector.zig");

const StridedArray = @import("stridedarray.zig").StridedArray;
const Vector       = @import("vector.zig"      ).Vector      ;

const asfloat = @import("helper.zig").asfloat;
const bitrev  = @import("helper.zig").bitrev ;
const prod    = @import("helper.zig").prod   ;

pub fn fft(comptime T: type, out: StridedArray(Complex(T)), in: StridedArray(Complex(T)), factor: i32) void {
    const n = in.len; const logn: u6 = @intCast(std.math.log2(n)); for (0..n) |i| out.ptr(bitrev(i, logn)).* = in.at(i);

    for (0..logn) |i| {

        const m = std.math.pow(usize, 2, i + 1); const mix = std.math.complex.exp(Complex(T).init(0, asfloat(T, factor) * 2 * std.math.pi / asfloat(T, m)));

        for (0..n / m) |j| {

            var omega = Complex(T).init(1, 0);

            for (0..m / 2) |k| {

                const t = omega.mul(out.at(j * m + k + m / 2)); const u = out.at(j * m + k);

                out.ptr(j * m + k        ).* = u.add(t);
                out.ptr(j * m + k + m / 2).* = u.sub(t);

                omega = omega.mul(mix);
            }
        }
    }

    if (factor > 0) for (0..out.len) |i| {
        out.ptr(i).* = out.at(i).div(Complex(T).init(asfloat(T, out.len), 0));
    };
}

pub fn fftn(comptime T: type, out: []Complex(T), in: []Complex(T), shape: []const usize, factor: i32) void {
    @memcpy(out, in); const sprod = prod(usize, shape); var stride: usize = 1;

    for (0..shape.len) |i| {

        @memcpy(in, out);

        for (0..sprod / shape[shape.len - i - 1]) |j| {

            const in_sa  = StridedArray(Complex(T)){.data = in,  .len = shape[shape.len - i - 1], .stride = stride, .zero = j * sprod / stride / shape[shape.len - i - 1]};
            const out_sa = StridedArray(Complex(T)){.data = out, .len = shape[shape.len - i - 1], .stride = stride, .zero = j * sprod / stride / shape[shape.len - i - 1]};

            fft(T, out_sa, in_sa, factor);
        }

        stride *= shape[shape.len - i - 1];
    }
}
