const std = @import("std"); const Complex = std.math.Complex; const gsl_fft = @cImport(@cInclude("gsl/gsl_fft_complex.h"));

const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;

pub fn fft(comptime T: type, array: []Complex(T), factor: i32, GSLFT: *gsl_fft.gsl_fft_complex_wavetable, GSLFW: *gsl_fft.gsl_fft_complex_workspace) void {
    if (factor == -1) _ = gsl_fft.gsl_fft_complex_forward (&array[0].re, 1, array.len, GSLFT, GSLFW);
    if (factor ==  1) _ = gsl_fft.gsl_fft_complex_backward(&array[0].re, 1, array.len, GSLFT, GSLFW);

    if (factor > 0) for (0..array.len) |i| {array[i] = array[i].div(Complex(T).init(asfloat(T, array.len), 0));};
}

test "fft" {
    var A = try Vector(Complex(f64)).init(4, std.testing.allocator); defer A.deinit();

    A.ptr(0).* = Complex(f64).init(3, 0);
    A.ptr(1).* = Complex(f64).init(6, 0);
    A.ptr(2).* = Complex(f64).init(1, 0);
    A.ptr(3).* = Complex(f64).init(8, 0);

    const GSLFT = gsl_fft.gsl_fft_complex_wavetable_alloc(4);
    const GSLFW = gsl_fft.gsl_fft_complex_workspace_alloc(4);

    fft(f64, A.data, -1, GSLFT, GSLFW);

    for (0..A.data.len) |i| {
        std.debug.print("{d}: {d:12.6} {d:12.6}\n", .{i + 1, A.data[i].re, A.data[i].im});
    }

    try std.testing.expect(true);
}
