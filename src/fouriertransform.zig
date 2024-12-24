const std = @import("std"); const Complex = std.math.Complex;

const vec = @import("vector.zig");

const StridedArray = @import("stridedarray.zig").StridedArray;
const Vector       = @import("vector.zig"      ).Vector      ;

const asfloat = @import("helper.zig").asfloat;
const bitrev  = @import("helper.zig").bitrev ;
const prod    = @import("helper.zig").prod   ;

pub fn fft(comptime T: type, arr: StridedArray(Complex(T)), factor: i32) !void {
    const n = arr.len; const logn: u6 = @intCast(std.math.log2(n));

    if (std.math.pow(usize, 2, @intCast(logn)) != n) {
        return error.InvalidFourierTransformLength;
    }

    for (0..n) |i| {

        const j = bitrev(i, logn);

        if (i < j) {
            std.mem.swap(Complex(T), arr.ptr(j), arr.ptr(i));
        }
    }

    for (0..logn) |i| {

        const m = std.math.pow(usize, 2, i + 1); const mix = std.math.complex.exp(Complex(T).init(0, asfloat(T, factor) * 2 * std.math.pi / asfloat(T, m)));

        for (0..n / m) |j| {

            var omega = Complex(T).init(1, 0);

            for (0..m / 2) |k| {

                const t = omega.mul(arr.at(j * m + k + m / 2)); const u = arr.at(j * m + k);

                arr.ptr(j * m + k        ).* = u.add(t);
                arr.ptr(j * m + k + m / 2).* = u.sub(t);

                omega = omega.mul(mix);
            }
        }
    }

    if (factor > 0) for (0..n) |i| {
        arr.ptr(i).* = arr.at(i).div(Complex(T).init(asfloat(T, n), 0));
    };
}

pub fn fftn(comptime T: type, arr: []Complex(T), shape: []const usize, factor: i32) !void {
    const sprod = prod(usize, shape); var stride: usize = 1;

    if (std.math.pow(usize, shape[0], shape.len) != sprod) {
        return error.InvalidFourierTransformShape;
    }

    for (0..shape.len) |i| {

        for (0..sprod / shape[shape.len - i - 1]) |j| {

            var offset: usize = 0; var index: usize = 0;

            for (0..shape.len) |k| if (k != i) {
                offset += (j / std.math.pow(usize, shape[k], shape.len - index - 2) % shape[k]) * std.math.pow(usize, shape[k], k); index += 1;
            };

            try fft(T, StridedArray(Complex(T)){.data = arr, .len = shape[shape.len - i - 1], .stride = stride, .zero = offset}, factor);
        }

        stride *= shape[shape.len - i - 1];
    }
}
