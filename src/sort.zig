//! Module for sorting files

const std = @import("std");

const hlp = @import("helper.zig");
const inp = @import("input.zig" );
const mat = @import("matrix.zig");

const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;

/// The main function for sorting files
pub fn run(comptime T: type, opt: inp.SortOptions(T), print: bool, allocator: std.mem.Allocator) !void {
    var timer = try std.time.Timer.start(); var A = try mat.read(T, opt.input, allocator);

    if (print) try hlp.print(std.fs.File.stdout(), "\nREADING THE ARRAY OF {d} NUMBERS TOOK {D}", .{A.data.len, timer.read()});

    timer = try std.time.Timer.start();

    if (std.mem.eql(u8, opt.algorithm, "bubble")) {bubbleMatrix(T, &A, opt.column);}

    else return error.InvalidSortingAlgorithm;

    if (print and opt.print) {try hlp.print(std.fs.File.stdout(), "\n\nSORTED MATRIX\n", .{}); try A.print(std.fs.File.stdout());}

    if (print) try hlp.print(std.fs.File.stdout(), "\nSORTING THE ARRAY OF {d} NUMBERS TOOK {D}\n", .{A.data.len, timer.read()});

    timer = try std.time.Timer.start();

    if (opt.output != null) try A.write(opt.output.?);

    if (print and opt.output != null) try hlp.print(std.fs.File.stdout(), "WRITING THE ARRAY OF {d} NUMBERS TOOK {D}\n", .{A.data.len, timer.read()});
}

// The bubble sort algorithm for matrices.
pub fn bubbleMatrix(comptime T: type, data: *Matrix(T), column: usize) void {
    for (0..data.rows) |_| {
        for (1..data.rows) |i| {
            if (data.at(i, column) < data.at(i - 1, column)) {
                for (0..data.cols) |j| {
                    std.mem.swap(T, data.ptr(i, j), data.ptr(i - 1, j));
                }
            }
        }
    }
}
