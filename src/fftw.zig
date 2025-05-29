//! This file contains wrappers for the FFTW3 library.

const std = @import("std"); const Complex = std.math.Complex;

const fftw = @cImport(@cInclude("fftw3.h"));

pub fn fftwnd(in: []Complex(f64), shape: []const u32, factor: i32) void {
    const plan = fftw.fftw_plan_dft(@intCast(shape.len), &@as([]const i32, @ptrCast(shape))[0], &@as([][2]f64, @ptrCast(in))[0], &@as([][2]f64, @ptrCast(in))[0], factor, fftw.FFTW_ESTIMATE);

    fftw.fftw_execute(plan); fftw.fftw_destroy_plan(plan);

    if (factor > 0) for (0..in.len) |i| {
        in[i] = in[i].div(Complex(f64).init(@as(f64, @floatFromInt(in.len)), 0));
    };
}
