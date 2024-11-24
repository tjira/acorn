const std = @import("std"); const Complex = std.math.Complex; const gsl_fft = @cImport(@cInclude("gsl/gsl_fft_complex.h"));

const asfloat = @import("helper.zig").asfloat;

pub fn fft(comptime T: type, array: []Complex(T), factor: i32, GSLFT: *gsl_fft.gsl_fft_complex_wavetable, GSLFW: *gsl_fft.gsl_fft_complex_workspace) void {
    if (factor == -1) _ = gsl_fft.gsl_fft_complex_forward (&array[0].re, 1, array.len, GSLFT, GSLFW);
    if (factor ==  1) _ = gsl_fft.gsl_fft_complex_backward(&array[0].re, 1, array.len, GSLFT, GSLFW);

    if (factor > 0) for (0..array.len) |i| {array[i] = array[i].div(Complex(T).init(asfloat(T, array.len), 0));};
}
