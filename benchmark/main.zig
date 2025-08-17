const std = @import("std"); const builtin = @import("builtin"); const acorn = @import("acorn");

const contract = @import("contract.zig");
const dgees    = @import("dgees.zig"   );
const dgemm    = @import("dgemm.zig"   );
const dsyevd   = @import("dsyevd.zig"  );

/// Function to benchmark a given object with a specific title.
pub fn benchmark(title: []const u8, object: anytype, D: usize, N: usize, print: bool, allocator: std.mem.Allocator) !void {

    const digits = std.math.log10(N) + 1;

    var timer = try std.time.Timer.start();

    var args = try object.init(allocator, D); defer args.deinit();

    try std.io.getStdOut().writer().print("{s} BENCHMARK (SIZE: {d}, SAMPLES: {d}, INIT: {s})\n\n", .{title, D, N, std.fmt.fmtDuration(timer.read())});

    var times = try allocator.alloc(f64, N); defer allocator.free(times);

    for (0..N) |i| {

        timer.reset();

        try std.io.getStdOut().writer().print("{[i]d:[width]}/{[n]d}) RANDOMIZATION: ", .{.i = i + 1, .width = digits, .n = N});

        try args.randomize(i);

        try std.io.getStdOut().writer().print("{s:9}, EXECUTION: ", .{std.fmt.fmtDuration(timer.read())  });

        timer.reset();

        try args.function();

        times[i] = acorn.helper.asfloat(f64, timer.read());

        try std.io.getStdOut().writer().print("{s:9}\n", .{std.fmt.fmtDuration(@as(u64, @intFromFloat(times[i])))});

        if (print) try args.print();
    }

    const time_mean   = @as(u64, @intFromFloat(acorn.math.mean(  f64, times) catch 0));
    const time_stddev = @as(u64, @intFromFloat(acorn.math.stddev(f64, times) catch 0));

    try std.io.getStdOut().writer().print("\nAVERAGE {s} TIME: {s} ± {d}\n", .{title, std.fmt.fmtDuration(time_mean), std.fmt.fmtDuration(time_stddev)});
}

/// Main function for the benchmarks
pub fn main() !void {
    if (std.posix.getenv("OMP_NUM_THREADS") == null) {
        return error.UndefinedNumberOfThreads;
    }

    var gpa = std.heap.GeneralPurposeAllocator(.{}){}; const allocator = gpa.allocator();

    defer {
        if (gpa.deinit() == .leak) std.debug.panic("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    var name = std.ArrayList(u8).init(allocator); var D: usize = 2; var N: usize = 1; defer name.deinit();

    try parse(&name, &D, &N, allocator);

    if      (std.mem.eql(u8, name.items, "contract")) {try benchmark("TENSOR CONTRACTION",               contract.Args, D, N, false, allocator);}
    else if (std.mem.eql(u8, name.items, "dgees"   )) {try benchmark("REAL SCHUR DECOMPOSITION",         dgees.Args,    D, N, false, allocator);}
    else if (std.mem.eql(u8, name.items, "dgemm"   )) {try benchmark("MATRIX MULTIPLICATION",            dgemm.Args,    D, N, false, allocator);}
    else if (std.mem.eql(u8, name.items, "dsyevd"   )) {try benchmark("SYMMETRIC MATRIX DIAGONALIZATION", dsyevd.Args,   D, N, false, allocator);}

    else return error.InvalidBenchmarkName;
}

/// Parser for the arguments, where first argument is difficulty and second is number of samples.
pub fn parse(name: *std.ArrayList(u8), D: *usize, N: *usize, allocator: std.mem.Allocator) !void {
    var argv = try std.process.argsWithAllocator(allocator); defer argv.deinit(); var argc: usize = 0;

    _ = argv.next(); while (argv.next()) |arg| {

        if (argc == 0) {
            try name.appendSlice(arg);
        }

        else if (argc == 1) {
            D.* = try std.fmt.parseInt(usize, arg, 10);
        }

        else if (argc == 2) {
            N.* = try std.fmt.parseInt(usize, arg, 10);
        }

        else {
            return error.TooManyArguments;
        }

        argc += 1;
    }

    if (argc == 0) return error.NoBenchmarkNameProvided;

    if (N.* < 1) {return error.InvalidNumberOfSamples;} if (D.* < 1) {return error.InvalidDifficulty;}
}
