const std = @import("std"); const acorn = @import("acorn");

const benchmark = @import("benchmark.zig");

const Args = struct {
    A: acorn.Matrix(f64), B: acorn.Matrix(f64), C: acorn.Matrix(f64),

    pub fn init(allocator: std.mem.Allocator, D: usize) !Args {
        return .{
            .A = try acorn.Matrix(f64).init(D, D, allocator),
            .B = try acorn.Matrix(f64).init(D, D, allocator),
            .C = try acorn.Matrix(f64).init(D, D, allocator),
        };
    }

    pub fn deinit(self: *Args) void {
        self.A.deinit(); self.B.deinit(); self.C.deinit();
    }

    pub fn function(self: *Args) !void {
        acorn.cwrapper.Blas(f64).dgemm(&self.C, self.A, false, self.B, false);
    }

    pub fn print(self: Args) !void {
        std.debug.print("MATRIX A:\n", .{}); try self.A.print(std.io.getStdOut().writer());
        std.debug.print("MATRIX B:\n", .{}); try self.B.print(std.io.getStdOut().writer());
        std.debug.print("MATRIX C:\n", .{}); try self.C.print(std.io.getStdOut().writer());
    }

    pub fn randomize(self: *Args, seed: usize) !void {
        self.A.randn(0, 1, 2 * seed + 0);
        self.B.randn(0, 1, 2 * seed + 1);
    }
};

pub fn main() !void {
    var D: usize = 2; try benchmark.benchmark("MATRIX MULTIPLICATION", Args, &D, false);
}
