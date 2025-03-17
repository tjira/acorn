//! Module for sorting files

const std = @import("std");

const inp = @import("input.zig" );
const mat = @import("matrix.zig");

const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;

/// The main function for sorting files
pub fn run(comptime T: type, opt: inp.SortOptions(), print: bool, allocator: std.mem.Allocator) !void {
    var A = try mat.read(T, opt.input, allocator); var timer = try std.time.Timer.start();

    if (std.mem.eql(u8, opt.algorithm, "bubble")) {bubble(T, A.data);}

    else return error.InvalidSortingAlgorithm;

    if (print) try std.io.getStdOut().writer().print("\nSORTING THE ARRAY OF {d} NUMBERS TOOK {}\n", .{A.data.len, std.fmt.fmtDuration(timer.read())});

    try A.write(opt.output);
}

/// The bubble sort algorithm
pub fn bubble(comptime T: type, data: []T) void {
    var sorted = true;

    for (0..data.len) |_| {

        for (1..data.len) |i| {
            if (data[i] < data[i - 1]) {std.mem.swap(T, &data[i], &data[i - 1]); sorted = false;}
        }

        if (sorted) return;
    }
}
