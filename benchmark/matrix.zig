const std = @import("std");

const helper = @import("acorn").helper;
const matrix = @import("acorn").matrix;

const Matrix = @import("acorn").Matrix;

const allocator = std.heap.page_allocator;

pub fn matrix_mm(samples: u64) ![2]f64 {
    var A = try matrix.Matrix(f64).init(700, 700, allocator); defer A.deinit();
    var B = try matrix.Matrix(f64).init(700, 700, allocator); defer B.deinit();
    var C = try matrix.Matrix(f64).init(700, 700, allocator); defer C.deinit();

    A.randn(0, 1, 0); B.randn(0, 1, 1);

    const times = try allocator.alloc(f64, samples); defer allocator.free(times);

    for (0..samples) |i| {
        var timer = try std.time.Timer.start(); matrix.mm(f64, &C, A, B); times[i] = helper.asfloat(f64, timer.read()) / 1e6;
    }

    return .{helper.mean(f64, times), helper.stdev(f64, times)};
}

pub fn matrix_eigh(samples: u64) ![2]f64 {
    var A  = try matrix.Matrix(f64).init(40, 40, allocator); defer  A.deinit(); 
    var UA = try matrix.Matrix(f64).init(40, 40, allocator); defer UA.deinit();
    var UC = try matrix.Matrix(f64).init(40, 40, allocator); defer UC.deinit();
    var T1 = try matrix.Matrix(f64).init(40, 40, allocator); defer T1.deinit();
    var T2 = try matrix.Matrix(f64).init(40, 40, allocator); defer T2.deinit();

    A.randn(0, 1, 0); matrix.transpose(f64, &T1, A); matrix.add(f64, &A, A, T1);

    const times = try allocator.alloc(f64, samples); defer allocator.free(times);

    for (0..samples) |i| {
        var timer = try std.time.Timer.start(); matrix.eigh(f64, &UA, &UC, A, &T1, &T2); times[i] = helper.asfloat(f64, timer.read()) / 1e6;
    }

    return .{helper.mean(f64, times), helper.stdev(f64, times)};
}
