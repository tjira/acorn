const std = @import("std"); const builtin = @import("builtin"); const acorn = @import("acorn");

const contract = @import("contract.zig");
const dgees    = @import("dgees.zig"   );
const dgemm    = @import("dgemm.zig"   );
const dgesdd   = @import("dgesdd.zig"  );
const dsyevd   = @import("dsyevd.zig"  );

/// Function to benchmark a given object with a specific title.
pub fn benchmark(title: []const u8, object: anytype, D: usize, N: usize, S: usize, print: bool, allocator: std.mem.Allocator) !void {

    const digits = std.math.log10(N) + 1;

    var timer = try std.time.Timer.start();

    var args = try object.init(allocator, D); defer args.deinit();

    try acorn.helper.print(std.fs.File.stdout(), "{s} BENCHMARK (SIZE: {d}, SAMPLES: {d}, INIT: {D})\n\n", .{title, D, N, timer.read()});

    var times = try allocator.alloc(f64, N); defer allocator.free(times);

    for (0..N) |i| {

        timer.reset();

        try acorn.helper.print(std.fs.File.stdout(), "{[i]d:[width]}/{[n]d}) RANDOMIZATION: ", .{.i = i + 1, .width = digits, .n = N});

        try args.randomize(i + S);

        try acorn.helper.print(std.fs.File.stdout(), "{D:9}, EXECUTION: ", .{timer.read()});

        timer.reset();

        try args.function();

        times[i] = acorn.helper.asfloat(f64, timer.read());

        try acorn.helper.print(std.fs.File.stdout(), "{D:9}\n", .{@as(u64, @intFromFloat(times[i]))});

        if (print) try args.print();
    }

    const time_mean   = @as(u64, @intFromFloat(acorn.math.mean(  f64, times) catch 0));
    const time_stddev = @as(u64, @intFromFloat(acorn.math.stddev(f64, times) catch 0));

    try acorn.helper.print(std.fs.File.stdout(), "\nAVERAGE {s} TIME: {D} ± {D}\n", .{title, time_mean, time_stddev});
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

    var name = std.ArrayList(u8){}; var D: usize = 2; var N: usize = 1; var S: usize = 0; var P: u1 = 0; defer name.deinit(allocator);

    try parse(&name, &D, &N, &S, &P, allocator);

    if      (std.mem.eql(u8, name.items, "contract")) {try benchmark("TENSOR CONTRACTION",               contract.Args, D, N, S, P == 1, allocator);}
    else if (std.mem.eql(u8, name.items, "dgees"   )) {try benchmark("REAL SCHUR DECOMPOSITION",         dgees.Args,    D, N, S, P == 1, allocator);}
    else if (std.mem.eql(u8, name.items, "dgemm"   )) {try benchmark("MATRIX MULTIPLICATION",            dgemm.Args,    D, N, S, P == 1, allocator);}
    else if (std.mem.eql(u8, name.items, "dgesdd"  )) {try benchmark("SVD DECOMPOSITION",                dgesdd.Args,   D, N, S, P == 1, allocator);}
    else if (std.mem.eql(u8, name.items, "dsyevd"  )) {try benchmark("SYMMETRIC MATRIX DIAGONALIZATION", dsyevd.Args,   D, N, S, P == 1, allocator);}

    else return error.InvalidBenchmarkName;
}

/// Parser for the arguments, where first argument is difficulty and second is number of samples.
pub fn parse(name: *std.ArrayList(u8), D: *usize, N: *usize, S: *usize, P: *u1, allocator: std.mem.Allocator) !void {
    var argv = try std.process.argsWithAllocator(allocator); defer argv.deinit(); var argc: usize = 0;

    _ = argv.next(); while (argv.next()) |arg| {

        if (argc == 0) {
            try name.appendSlice(allocator, arg);
        }

        else if (argc == 1) {
            D.* = try std.fmt.parseInt(usize, arg, 10);
        }

        else if (argc == 2) {
            N.* = try std.fmt.parseInt(usize, arg, 10);
        }

        else if (argc == 3) {
            S.* = try std.fmt.parseInt(usize, arg, 10);
        }

        else if (argc == 4) {
            P.* = try std.fmt.parseInt(u1, arg, 10);
        }

        else {
            return error.TooManyArguments;
        }

        argc += 1;
    }

    if (argc == 0) return error.NoBenchmarkNameProvided;

    if (N.* < 1) {return error.InvalidNumberOfSamples;} if (D.* < 1) {return error.InvalidDifficulty;}
}
