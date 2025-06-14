const std = @import("std"); const builtin = @import("builtin"); const acorn = @import("acorn");

pub fn main() !void {
    if (std.posix.getenv("OMP_NUM_THREADS") == null) {
        return error.UndefinedNumberOfThreads;
    }

    std.debug.print("ZIG VERSION: {}, THREADS: {s}\n", .{builtin.zig_version, std.posix.getenv("OMP_NUM_THREADS").?});

    var gpa = std.heap.GeneralPurposeAllocator(.{}){}; const allocator = gpa.allocator();

    defer {
        if (gpa.deinit() == .leak) std.debug.panic("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    var argv = try std.process.argsWithAllocator(allocator); defer argv.deinit();

    var argc: usize = 0; var N: usize = 2; var R: usize = 4;

    _ = argv.next(); while (argv.next()) |arg| {

        if (argc == 0) {
            R = try std.fmt.parseInt(usize, arg, 10);
        }

        else if (argc == 1) {
            N = try std.fmt.parseInt(usize, arg, 10);
        }

        else {
            return error.TooManyArguments;
        }

        argc += 1;
    }

    if (N < 2) {return error.InvalidNumberOfMatrices;} if (R < 1) {return error.InvalidMatrixRank;}

    var micro_timer = try std.time.Timer.start();

    var A  = try acorn.Matrix(f64).init(R, R, allocator); defer  A.deinit();
    var AC = try acorn.Matrix(f64).init(R, R, allocator); defer AC.deinit();
    var AJ = try acorn.Matrix(f64).init(R, R, allocator); defer AJ.deinit();

    var results = try allocator.alloc(f64, N); defer allocator.free(results);

    std.debug.print("\n3 RANK {d} MATRICES ALLOCATION TIME: {s}\n", .{R, std.fmt.fmtDuration(micro_timer.read())});

    var macro_timer = try std.time.Timer.start();

    for (0..N) |i| {

        micro_timer = try std.time.Timer.start();

        A.randn(0, 1, i); try A.symmetrize();

        const time_gen = std.fmt.fmtDuration(micro_timer.read());

        micro_timer = try std.time.Timer.start();

        try acorn.cwrapper.dsyevd(&AJ, &AC, A);

        results[i] = acorn.helper.asfloat(f64, micro_timer.read());

        const time_dsyevd = std.fmt.fmtDuration(micro_timer.read());

        std.debug.print("{s}{d}/{d}: MATRIX GENERATION AND DIAGONALIZATION TIMING: {s}/{s}\n", .{if (i == 0) "\n" else "", i + 1, N, time_gen, time_dsyevd});
    }

    const result_mean   = @as(u64, @intFromFloat(try acorn.math.mean(  f64, results)));
    const result_stddev = @as(u64, @intFromFloat(try acorn.math.stddev(f64, results)));

    std.debug.print("\nESTIMATED DIAGONALIZATION TIME FOR RANK {d} MATRICES: {s}±{d}\n", .{R, std.fmt.fmtDuration(result_mean), std.fmt.fmtDuration(result_stddev)});

    std.debug.print("\nTOTAL EXECUTION TIME: {s}\n", .{std.fmt.fmtDuration(macro_timer.read())});
}
