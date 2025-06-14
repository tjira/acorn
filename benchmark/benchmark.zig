const std = @import("std"); const builtin = @import("builtin"); const acorn = @import("acorn");

/// Function to benchmark a given object with a specific title.
pub fn benchmark(title: []const u8, object: anytype, D: *usize, print: bool) !void {
    if (std.posix.getenv("OMP_NUM_THREADS") == null) {
        return error.UndefinedNumberOfThreads;
    }

    std.debug.print("ZIG VERSION: {}, THREADS: {s}\n\n", .{builtin.zig_version, std.posix.getenv("OMP_NUM_THREADS").?});

    var gpa = std.heap.GeneralPurposeAllocator(.{}){}; const allocator = gpa.allocator();

    defer {
        if (gpa.deinit() == .leak) std.debug.panic("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    var N: usize = 1; try parse(D, &N, allocator); const digits = std.math.log10(N) + 1;

    var args = try object.init(allocator, D.*); defer args.deinit();

    var times = try allocator.alloc(f64, N); defer allocator.free(times);

    for (0..N) |i| {

        try args.randomize(i);

        std.debug.print("{[i]d:[width]}/{[n]d}: {[title]s} TIME: ", .{.i = i + 1, .width = digits, .n = N, .title = title});

        var timer = try std.time.Timer.start();

        try args.function();

        times[i] = acorn.helper.asfloat(f64, timer.read());

        std.debug.print("{s}\n", .{std.fmt.fmtDuration(@as(u64, @intFromFloat(times[i])))});

        if (print) try args.print();
    }

    const time_mean   = @as(u64, @intFromFloat(acorn.math.mean(  f64, times) catch 0));
    const time_stddev = @as(u64, @intFromFloat(acorn.math.stddev(f64, times) catch 0));

    std.debug.print("\nAVERAGE {s} TIME: {s} ± {d}\n", .{title, std.fmt.fmtDuration(time_mean), std.fmt.fmtDuration(time_stddev)});
}

/// Parser for the arguments, where first argument is difficulty and second is number of samples.
pub fn parse(D: *usize, N: *usize, allocator: std.mem.Allocator) !void {
    var argv = try std.process.argsWithAllocator(allocator); defer argv.deinit(); var argc: usize = 0;

    _ = argv.next(); while (argv.next()) |arg| {

        if (argc == 0) {
            D.* = try std.fmt.parseInt(usize, arg, 10);
        }

        else if (argc == 1) {
            N.* = try std.fmt.parseInt(usize, arg, 10);
        }

        else {
            return error.TooManyArguments;
        }

        argc += 1;
    }

    if (N.* < 1) {return error.InvalidNumberOfSamples;} if (D.* < 1) {return error.InvalidDifficulty;}
}
