//! Fourier transform module.

const std = @import("std");

const mth = @import("math.zig"  );
const vec = @import("vector.zig");

const StridedArray = @import("strided_array.zig").StridedArray;
const Vector       = @import("vector.zig"       ).Vector      ;

const asfloat = @import("helper.zig").asfloat;
const bitrev  = @import("helper.zig").bitrev ;

/// Fast Fourier transform for a one-dimensional array. The factor argument is the value in the exponent of the Fourier transform. Factor -1 corresponds to the forward Fourier transform, while factor 1 corresponds to the inverse Fourier transform.
pub fn fft(comptime T: type, arr: StridedArray(std.math.Complex(T)), factor: i32) !void {
    const n = arr.len; const logn: u6 = @intCast(std.math.log2(n));

    if (std.math.pow(usize, 2, @intCast(logn)) != n) {
        return error.InvalidFourierTransformLength;
    }

    for (0..n) |i| {

        const j = bitrev(i, logn);

        if (i < j) {
            std.mem.swap(std.math.Complex(T), arr.ptr(j), arr.ptr(i));
        }
    }

    for (0..logn) |i| {

        const m = std.math.pow(usize, 2, i + 1); const mix = std.math.complex.exp(std.math.Complex(T).init(0, asfloat(T, factor) * 2 * std.math.pi / asfloat(T, m)));

        for (0..n / m) |j| {

            var omega = std.math.Complex(T).init(1, 0);

            for (0..m / 2) |k| {

                const t = omega.mul(arr.at(j * m + k + m / 2)); const u = arr.at(j * m + k);

                arr.ptr(j * m + k        ).* = u.add(t);
                arr.ptr(j * m + k + m / 2).* = u.sub(t);

                omega = omega.mul(mix);
            }
        }
    }

    if (factor > 0) for (0..n) |i| {
        arr.ptr(i).* = arr.at(i).div(std.math.Complex(T).init(asfloat(T, n), 0));
    };
}

/// Fast Fourier transform for an n-dimensional array. The factor argument is the value in the exponent of the Fourier transform. Factor -1 corresponds to the forward Fourier transform, while factor 1 corresponds to the inverse Fourier transform.
pub fn fftn(comptime T: type, arr: []std.math.Complex(T), shape: []const usize, factor: i32) !void {
    const N = mth.prod(usize, shape);

    if (std.math.pow(usize, shape[0], shape.len) != N) {
        return error.InvalidFourierTransformShape;
    }

    for (0..shape.len) |i| {

        const stride = mth.prod(usize, shape[0..i]);

        for (0..N / shape[i]) |j| {

            var offset: usize = 0; var index: usize = 0;

            for (0..shape.len) |k| if (k != i) {

                const block = mth.prod(usize, shape[k + 1..]) / (if (k < i) shape[i] else 1);

                offset += (j / block % shape[k]) * mth.prod(usize, shape[0..k]); index += 1;
            };

            try fft(T, StridedArray(std.math.Complex(T)){.data = arr, .len = shape[i], .stride = stride, .zero = offset}, factor);
        }
    }
}
