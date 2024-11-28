const std = @import("std"); const Complex = std.math.Complex; const gsl_fft = @cImport(@cInclude("gsl/gsl_fft_complex.h"));

const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;

pub fn fft(comptime T: type, array: []Complex(T), factor: i32, GSLFT: *gsl_fft.gsl_fft_complex_wavetable, GSLFW: *gsl_fft.gsl_fft_complex_workspace) void {
    if (factor == -1) _ = gsl_fft.gsl_fft_complex_forward (&array[0].re, 1, array.len, GSLFT, GSLFW);
    if (factor ==  1) _ = gsl_fft.gsl_fft_complex_backward(&array[0].re, 1, array.len, GSLFT, GSLFW);

    if (factor > 0) for (0..array.len) |i| {array[i] = array[i].div(Complex(T).init(asfloat(T, array.len), 0));};
}

pub fn myfft(comptime T: type, out: []Complex(T), in: []Complex(T), factor: i32) !void {
    // const n = in.len;

    if (in.len == 1) {out[0] = in[0]; return;}

    var even = try std.heap.page_allocator.alloc(Complex(T), in.len / 2); defer std.heap.page_allocator.free(even);
    var odd  = try std.heap.page_allocator.alloc(Complex(T), in.len / 2); defer std.heap.page_allocator.free(odd );

    for (0..in.len / 2) |i| {
        even[i] = in[2 * i + 0];
        odd[i]  = in[2 * i + 1];
    }

    const fft_even = try std.heap.page_allocator.alloc(Complex(T), in.len / 2); defer std.heap.page_allocator.free(fft_even);
    const fft_odd  = try std.heap.page_allocator.alloc(Complex(T), in.len / 2); defer std.heap.page_allocator.free(fft_odd );

    try myfft(T, fft_even, even, factor);
    try myfft(T, fft_odd,  odd,  factor);

    for (0..in.len / 2) |k| {
        out[k             ] = fft_even[k].add(std.math.complex.exp(Complex(T).init(0, asfloat(T, factor) * 2 * std.math.pi * asfloat(T, k) / asfloat(T, in.len))).mul(fft_odd[k]));
        out[k + in.len / 2] = fft_even[k].sub(std.math.complex.exp(Complex(T).init(0, asfloat(T, factor) * 2 * std.math.pi * asfloat(T, k) / asfloat(T, in.len))).mul(fft_odd[k]));
    }
}

test "fft" {
    var A = try Vector(Complex(f64)).init(4, std.testing.allocator); defer A.deinit();
    var B = try Vector(Complex(f64)).init(4, std.testing.allocator); defer B.deinit();

    A.ptr(0).* = Complex(f64).init(3, 0);
    A.ptr(1).* = Complex(f64).init(6, 0);
    A.ptr(2).* = Complex(f64).init(1, 0);
    A.ptr(3).* = Complex(f64).init(8, 0);

    // const GSLFT = gsl_fft.gsl_fft_complex_wavetable_alloc(4);
    // const GSLFW = gsl_fft.gsl_fft_complex_workspace_alloc(4);

    // fft(f64, A.data, -1, GSLFT, GSLFW);
    try myfft(f64, B.data, A.data, -1); @memcpy(A.data, B.data);

    for (0..A.data.len) |i| {
        std.debug.print("{d}: {d:12.6} {d:12.6}\n", .{i + 1, A.data[i].re, A.data[i].im});
    }

    try std.testing.expect(true);
}
