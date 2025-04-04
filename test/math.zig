const std = @import("std");

const math = @import("acorn").math;

const expect = @import("main.zig").expect;
const log    = @import("main.zig").log   ;

const allocator = std.testing.allocator;

test "math_boys" {
    try expect(math.boys(f64,  0.01,  0.1), 0.83021792561502);
    try expect(math.boys(f64,  0.10,  0.1), 0.80302216920873);
    try expect(math.boys(f64,  1.00,  0.1), 0.59771744503671);
    try expect(math.boys(f64, 10.00,  0.1), 0.18703190059280);

    try expect(math.boys(f64,  0.01,  1.0), 0.33134045770977);
    try expect(math.boys(f64,  0.10,  1.0), 0.31402947299816);
    try expect(math.boys(f64,  1.00,  1.0), 0.18947234582049);
    try expect(math.boys(f64, 10.00,  1.0), 0.01401009952884);

    try expect(math.boys(f64,  0.01, 10.0), 0.04718625885185);
    try expect(math.boys(f64,  0.10, 10.0), 0.04346518972410);
    try expect(math.boys(f64,  1.00, 10.0), 0.01917293609131);
    try expect(math.boys(f64, 10.00, 10.0), 0.00000857837826);
}

test "math_gamma" {
    try expect(math.gamma(f64,   0.01),     99.43258511915057);
    try expect(math.gamma(f64,   0.10),      9.51350769866873);
    try expect(math.gamma(f64,   1.00),      1.00000000000000);
    try expect(math.gamma(f64,  10.00), 362880.00000000150000);
}

test "math_gammainc" {
    try expect(math.gammainc(f64,  0.01,  0.1), 0.66262125995448);
    try expect(math.gammainc(f64,  0.10,  0.1), 0.82755175958585);
    try expect(math.gammainc(f64,  1.00,  0.1), 0.97587265627367);
    try expect(math.gammainc(f64, 10.00,  0.1), 0.99999944520143);

    try expect(math.gammainc(f64,  0.01,  1.0), 0.00995016625083);
    try expect(math.gammainc(f64,  0.10,  1.0), 0.09516258196404);
    try expect(math.gammainc(f64,  1.00,  1.0), 0.63212055882856);
    try expect(math.gammainc(f64, 10.00,  1.0), 0.99995460007024);

    try expect(math.gammainc(f64,  0.01, 10.0), 0.00000000000000);
    try expect(math.gammainc(f64,  0.10, 10.0), 0.00000000000000);
    try expect(math.gammainc(f64,  1.00, 10.0), 0.00000011142548);
    try expect(math.gammainc(f64, 10.00, 10.0), 0.54207028552815);
}
