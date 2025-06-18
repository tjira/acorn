const std = @import("std"); const acorn = @import("acorn");

const benchmark = @import("benchmark.zig");

const Args = struct {
    A: acorn.Matrix(f64), AT: acorn.Matrix(f64), AQ: acorn.Matrix(f64), AJR: acorn.Matrix(f64), AJI: acorn.Matrix(f64),

    pub fn init(allocator: std.mem.Allocator, D: usize) !Args {
        return .{
            .A   = try acorn.Matrix(f64).init(D, D, allocator),
            .AT  = try acorn.Matrix(f64).init(D, D, allocator),
            .AQ  = try acorn.Matrix(f64).init(D, D, allocator),
            .AJR = try acorn.Matrix(f64).init(D, D, allocator),
            .AJI = try acorn.Matrix(f64).init(D, D, allocator),
        };
    }

    pub fn deinit(self: *Args) void {
        self.A.deinit(); self.AT.deinit(); self.AQ.deinit(); self.AJR.deinit(); self.AJI.deinit();
    }

    pub fn function(self: *Args) !void {
        try acorn.cwrapper.dgees(&self.AT, &self.AQ, self.A, &self.AJR, &self.AJI);
    }

    pub fn print(self: Args) !void {
        std.debug.print("MATRIX A:\n",        .{}); try self.A.print( std.io.getStdOut().writer());
        std.debug.print("SCHUR FORM:\n",      .{}); try self.AT.print(std.io.getStdOut().writer());
        std.debug.print("ORTHOGONAL FORM:\n", .{}); try self.AQ.print(std.io.getStdOut().writer());
    }

    pub fn randomize(self: *Args, seed: usize) !void {
        self.A.randn(0, 1, seed);
    }
};

pub fn main() !void {
    var D: usize = 2; try benchmark.benchmark("REAL SCHUR DECOMPOSITION", Args, &D, false);
}
