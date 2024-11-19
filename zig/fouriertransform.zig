const std = @import("std");

const fftw = @cImport(@cInclude("fftw3.h"));

const asfloat = @import("helper.zig").asfloat;

pub fn dft(comptime T: type, out: []std.math.Complex(T), in: []const std.math.Complex(T), shape: []const usize, factor: T) void {
    const size = shape[0]; var exp = std.math.Complex(T).init(0, 0); for (0..size) |i| out[i] = std.math.Complex(T).init(0, 0);

    for (0..size) |i| for (0..size) |j| {
        exp.im = factor * 2 * std.math.pi * asfloat(T, i * j) / asfloat(T, size); out[i] = out[i].add(std.math.complex.exp(exp).mul(in[j]));
    };

    if (factor > 0) for (0..size) |i| {out[i] = out[i].div(std.math.Complex(T).init(asfloat(T, size), 0));};
}

pub fn fft_fftw(comptime T: type, out: []std.math.Complex(T), in: []const std.math.Complex(T), shape: []const usize, factor: i32) !void {
    const size = shape[0];

    var fftw_in  = try std.ArrayList(fftw.fftw_complex).initCapacity(std.heap.page_allocator, size); defer  fftw_in.deinit();
    var fftw_out = try std.ArrayList(fftw.fftw_complex).initCapacity(std.heap.page_allocator, size); defer fftw_out.deinit();

    for (0..size) |i| try fftw_in .append(.{ in[i].re,  in[i].im});
    for (0..size) |i| try fftw_out.append(.{out[i].re, out[i].im});

    const plan = fftw.fftw_plan_dft_1d(@intCast(size), fftw_in.items.ptr, fftw_out.items.ptr, factor, fftw.FFTW_ESTIMATE); defer fftw.fftw_destroy_plan(plan);

    fftw.fftw_execute(plan); for (0..size) |i| out[i] = std.math.Complex(T).init(fftw_out.items[i][0], fftw_out.items[i][1]);

    if (factor > 0) for (0..size) |i| {out[i] = out[i].div(std.math.Complex(T).init(asfloat(T, size), 0));};
}
